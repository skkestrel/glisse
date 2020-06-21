#pragma once
#include <string>
#include <memory>
#include <utility>
#include <cstdlib>
#include <istream>
#include <ostream>
#include <functional>
#include "util.h"
#include "types.h"

using size_t = std::size_t;

namespace sr
{
namespace data
{
	/**
	 * This class represents a set of particle states.
	 * There is no data for death time or particle flags,
	 * so this class is used to represent a snapshot that
	 * makes up a single entry in a simulation track.
	 */
	struct HostParticleSnapshot
	{
		/** The number of particles that this snapshot contains.  */
		size_t n;

		/**
		 * The number of alive particles. The alive particles always come before
		 * dead particles in the array in order to speed up CUDA calls.
		 * However, when running in CPU-only mode, the array maybe out of order
		 * after integrator->step_particles_timeblock() is called and before resync() is called.
		 */
		size_t n_alive;

		/**
		 * The cartesian position vectors of the particles.
		 * These can be either barycentric or heliocentric.
		 * These may also contain the orbital elements a, e, i
		 * when a HostParticleSnapshot is returned from reading
		 * particle tracks.
		 */
		Vf64_3 r;

		/**
		 * The cartesian velocity vectors of the particles.
		 * These can be either barycentric or heliocentric.
		 * These may also contain the orbital elements O, o, f
		 * when a HostParticleSnapshot is returned from reading
		 * particle tracks.
		 */
		Vf64_3 v;

		/**
		 * The IDs of the particles.
		 */
		Vu32 id;

		/**
		 * Implements the gather operation on all of the particle arrays.
		 * The gather operation reorders the particle data in the order provided
		 * by the `indices` array.
		 * See `sr::data::gather<T>` for details on the gather operation.
		 */
		void gather(const std::vector<size_t>& indices, size_t begin, size_t length);

		/**
		 * Resizes all of the particle arrays to the given length.
		 */
		void resize(size_t length);

		/**
		 * Sorts all of the particles and their corresponding arrays
		 * by their ID. This function returns an array of indices that can
		 * be passed to a gather operation to reorder other data
		 * that particles may have.
		 */
		std::unique_ptr<std::vector<size_t>> sort_by_id(size_t begin, size_t length);

		/**
		 * Copies particle data to the given `sr::data::HostParticleSnapshot`,
		 * but only the particles with indices (not IDs) that the `filter` argument contains.
		 */
		void filter(const std::vector<size_t>& filter, HostParticleSnapshot& out) const;

		/** Default ctor. */
		inline HostParticleSnapshot() : n(0), n_alive(0) { }
		/** Ctor with size argument. */
		inline HostParticleSnapshot(size_t n_) : n(n_), n_alive(n_), r(n_), v(n_), id(n_) { }
	};

	struct HostPlanetSnapshot
	{
		/** The number of planets that this snapshot contains. */
		size_t n;

		/** The number of alive planets that this snapshot contains. Currently is equivalent to `n`. */
		size_t n_alive;

		/**
		 * The cartesian position vectors of the planets.
		 * These can be either barycentric or heliocentric.
		 * These may also contain the orbital elements a, e, i
		 * when a HostPlanetSnapshot is returned from reading
		 * particle tracks.
		 */
		Vf64_3 r;

		/**
		 * The cartesian veloicty vectors of the planets.
		 * These can be either barycentric or heliocentric.
		 * These may also contain the orbital elements O, o, f
		 * when a HostPlanetSnapshot is returned from reading
		 * particle tracks.
		 */
		Vf64_3 v;

		/**
		 * The IDs of the planets.
		 */
		Vu32 id;

		/**
		 * The masses of the planets.
		 */
		Vf64 m;

		/** Default ctor. */
		inline HostPlanetSnapshot() : n(0), n_alive(0) { }
		/** Ctor with size argument. */
		inline HostPlanetSnapshot(size_t n_) : n(n_), n_alive(n_), r(n_), v(n_), id(n_), m(n_) { }
	};


	struct HostParticlePhaseSpace
	{
		/** Basic particle data. */
		HostParticleSnapshot base;

		/** Gets the particle position array. */
		inline Vf64_3& r() { return base.r; }
		inline const Vf64_3& r() const { return base.r; }

		/** Gets the particle velocity array. */
		inline Vf64_3& v() { return base.v; }
		inline const Vf64_3& v() const { return base.v; }

		/** Gets the particle ID array. */
		inline Vu32& id() { return base.id; }
		inline const Vu32& id() const { return base.id; }

		/** Gets the number of particles. */
		inline size_t& n() { return base.n; }
		inline const size_t& n() const { return base.n; }

		/** Gets the number of alive particles. */
		inline size_t& n_alive() { return base.n_alive; }
		inline const size_t& n_alive() const { return base.n_alive; }

		/**
		 * Gets the particle death flags.
		 * Guide:
		 * High byte:
		 *   Planet ID that killed the particle
		 * Low byte:
		 *   0x80 Particle absorbed by planet
		 *   0x04 Kepler didn't converge
		 *   0x02 Out of bounds
		 */
		inline Vu16& deathflags() { return _deathflags; }
		inline const Vu16& deathflags() const { return _deathflags; }

		/**
		 * Gets the particle death time.
		 * Particles that remain alive have their death time set to the current time.
		 */
		inline Vf32& deathtime() { return _deathtime; }
		inline const Vf32& deathtime() const { return _deathtime; }

		/**
		 * CPU block integration mode only.
		 * Gets the index within the timeblock of the death time. This allows
		 * encounter integration to continue from the correct place.
		 */
		inline Vu32& deathtime_index() { return _deathtime_index; }
		inline const Vu32& deathtime_index() const { return _deathtime_index; }

		/** Default ctor. */
		inline HostParticlePhaseSpace() { }

		/** Ctor with size argument */
		inline HostParticlePhaseSpace(size_t siz) : base(siz), _deathflags(siz), _deathtime(siz) { }

		/** Execute stable partition on alive particles, i.e., `deathflags & 0x00FE = 0` */
		std::unique_ptr<std::vector<size_t>> stable_partition_alive(size_t begin, size_t length);

		/** Execute stable partition on unflagged particles, i.e., `deathflags & 0x00FF = 0`. */
		std::unique_ptr<std::vector<size_t>> stable_partition_unflagged(size_t begin, size_t length);

		/**
		 * Implements the gather operation on all of the particle arrays.
		 * The gather operation reorders the particle data in the order provided
		 * by the indices array.
		 * See sr::data::gather for details on the gather operation.
		 */
		void gather(const std::vector<size_t>& indices, size_t begin, size_t length);

		/**
		 * Sorts particle data by IDs. See `HostParticleSnapshot::sort_by_id`.
		 */
		std::unique_ptr<std::vector<size_t>> sort_by_id(size_t begin, size_t length);

		/**
		 * Filters particle data. See `HostParticleSnapshot::filter`.
		 */
		void filter(const std::vector<size_t>& filter, HostParticlePhaseSpace& out) const;

	private:

		Vu16 _deathflags;
		Vf32 _deathtime;

		Vu32 _deathtime_index;
	};

	struct HostPlanetPhaseSpace
	{
		/**
		 * Basic planet data.
		 */
		HostPlanetSnapshot base;

		/**
		 * Gets the planet position array.
		 */
		inline Vf64_3& r() { return base.r; }
		inline const Vf64_3& r() const { return base.r; }

		/**
		 * Gets the planet velocity array.
		 */
		inline Vf64_3& v() { return base.v; }
		inline const Vf64_3& v() const { return base.v; }

		/**
		 * Gets the planet ID array.
		 */
		inline Vu32& id() { return base.id; }
		inline const Vu32& id() const { return base.id; }

		/**
		 * Gets the planet mass array.
		 */
		inline Vf64& m() { return base.m; }
		inline const Vf64& m() const { return base.m; }

		/**
		 * Gets the planet count.
		 */
		inline size_t& n() { return base.n; }
		inline const size_t& n() const { return base.n; }

		/**
		 * Gets the alive planet count in the current timeblock.
		 */
		inline size_t& n_alive() { return base.n_alive; }
		inline const size_t& n_alive() const { return base.n_alive; }

		/**
		 * Gets the alive planet count for the previous timeblock.
		 * This is used in cases where the number of planets has changed.
		 */
		inline size_t& n_alive_old() { return _n_alive_old; }
		inline const size_t& n_alive_old() const { return _n_alive_old; }

		/**
		 * Gets various log of positions of the planets.
		 */
		inline sr::util::LogQuartet<Vf64_3>& r_log() { return _r_log; }
		inline const sr::util::LogQuartet<Vf64_3>& r_log() const { return _r_log; }

		/**
		 * Gets various log of velocities of the planets.
		 */
		inline sr::util::LogQuartet<Vf64_3>& v_log() { return _v_log; }
		inline const sr::util::LogQuartet<Vf64_3>& v_log() const { return _v_log; }

		inline HostPlanetPhaseSpace() { }
		inline HostPlanetPhaseSpace(size_t siz, size_t tb_size) :
			base(siz), _n_alive_old(siz)
		{
			_r_log = sr::util::LogQuartet<Vf64_3>(((n() - 1) * tb_size));
			_v_log = sr::util::LogQuartet<Vf64_3>(((n() - 1) * tb_size));
		}

		/**
		 * Swap data for the previous and current timeblocks.
		 */
		inline void swap_logs()
		{
			std::swap(n_alive(), n_alive_old());

			r_log().swap_logs();
			v_log().swap_logs();
		}

		/**
		 * Get the index in a log for a given timestep and planet ID.
		 */
		template<bool old>
		inline size_t log_index_at(size_t timestep, size_t planet_id) const
		{
			return timestep * ((old ? n_alive() : _n_alive_old) - 1) + planet_id - 1;
		}

	private:
		size_t _n_alive_old;

		sr::util::LogQuartet<Vf64_3> _r_log;
		sr::util::LogQuartet<Vf64_3> _v_log;
	};

	/** This struct represents data stored on the CPU. */
	struct HostData
	{
		/** Particle data. */
		HostParticlePhaseSpace particles;
		/** Planet data. */
		HostPlanetPhaseSpace planets;
		/**
		 * Frozen planet data.
		 * Since the planet data is always ahead
		 * of the particle data due to pre-integration,
		 * when a dump is needed, planet data is taken
		 * from the snapshot instead.
		 */
		HostPlanetSnapshot planets_snapshot;

		inline HostData() { }
	};

	/** Contains all the integration options. */
	struct Configuration
	{
		uint32_t max_kep;
		double t_0, t_f, dt, big_g;
		uint32_t num_thread;
		uint32_t tbsize, print_every, dump_every, track_every, energy_every, max_particle;
		double wh_ce_r1, wh_ce_r2;
		uint32_t wh_ce_n1, wh_ce_n2;
		uint32_t split_track_file;

		uint32_t resync_every;

		bool write_bary_track;

		double cull_radius;

		bool readmomenta, writemomenta, trackbinary, readsplit, writesplit, dumpbinary, writebinary, readbinary;

		std::string icsin, plin, hybridin, hybridout;
		std::string outfolder;

		Configuration();

		/**
		 * Create an output configuration for the current configuration.
		 */
		Configuration output_config() const;

		/**
		 * Create a dummy configuration for using data loading and writing in the data namespace.
		 */
		inline static Configuration create_dummy()
		{
			Configuration config;
			config.write_bary_track = false;
			config.tbsize = 0;
			config.wh_ce_n1 = 0;
			config.wh_ce_n2 = 0;

			return config;
		}
	};

	template<typename T>
	inline T reverse_2byte(T in)
	{
		static_assert(sizeof(T) == 2, "");
		T ret;
		char *inb = reinterpret_cast<char*>(&in);
		char *retb = reinterpret_cast<char*>(&ret);

		// swap the bytes into a temporary buffer
		retb[0] = inb[1];
		retb[1] = inb[0];

		return ret;
	}

	template<typename T>
	inline T reverse_4byte(T in)
	{
		static_assert(sizeof(T) == 4, "");
		T ret;
		char *inb = reinterpret_cast<char*>(&in);
		char *retb = reinterpret_cast<char*>(&ret);

		// swap the bytes into a temporary buffer
		retb[0] = inb[3];
		retb[1] = inb[2];
		retb[2] = inb[1];
		retb[3] = inb[0];

		return ret;
	}

	template<typename T>
	inline T reverse_8byte(T in)
	{
		static_assert(sizeof(T) == 8, "");
		T ret;
		char *inb = reinterpret_cast<char*>(&in);
		char *retb = reinterpret_cast<char*>(&ret);

		// swap the bytes into a temporary buffer
		retb[0] = inb[7];
		retb[1] = inb[6];
		retb[2] = inb[5];
		retb[3] = inb[4];
		retb[4] = inb[3];
		retb[5] = inb[2];
		retb[6] = inb[1];
		retb[7] = inb[0];

		return ret;
	}

	template<typename T>
	inline T reverse_bytes(T in) { (void) in; return T::unimplemented; }

	template<>
	inline int16_t reverse_bytes<int16_t>(int16_t in) { return reverse_2byte(in); }
	template<>
	inline uint16_t reverse_bytes<uint16_t>(uint16_t in) { return reverse_2byte(in); }
	template<>
	inline int32_t reverse_bytes<int32_t>(int32_t in) { return reverse_4byte(in); }
	template<>
	inline uint32_t reverse_bytes<uint32_t>(uint32_t in) { return reverse_4byte(in); }
	template<>
	inline int64_t reverse_bytes<int64_t>(int64_t in) { return reverse_8byte(in); }
	template<>
	inline uint64_t reverse_bytes<uint64_t>(uint64_t in) { return reverse_8byte(in); }
	template<>
	inline double reverse_bytes<double>(double in) { return reverse_8byte(in); }
	template<>
	inline float reverse_bytes<float>(float in) { return reverse_4byte(in); }

	inline bool is_int_little_endian()
	{
		uint32_t i = 0x01;
		return reinterpret_cast<uint8_t*>(&i)[0] == 0x01;
	}

	inline bool is_float_little_endian()
	{
		float i = 2;
		return reinterpret_cast<uint8_t*>(&i)[3] == 0x40;
	}

	inline bool is_double_little_endian()
	{
		double i = 2;
		return reinterpret_cast<uint8_t*>(&i)[7] == 0x40;
	}

	template<typename T>
	inline bool is_little_endian() { return T::unimplemented; }

	template<>
	inline bool is_little_endian<int16_t>() { return is_int_little_endian(); }
	template<>
	inline bool is_little_endian<uint16_t>() { return is_int_little_endian(); }
	template<>
	inline bool is_little_endian<int32_t>() { return is_int_little_endian(); }
	template<>
	inline bool is_little_endian<uint32_t>() { return is_int_little_endian(); }
	template<>
	inline bool is_little_endian<float>() { return is_float_little_endian(); }
	template<>
	inline bool is_little_endian<int64_t>() { return is_int_little_endian(); }
	template<>
	inline bool is_little_endian<uint64_t>() { return is_int_little_endian(); }
	template<>
	inline bool is_little_endian<double>() { return is_double_little_endian(); }

	template<typename T>
	inline T to_little_endian(T in)
	{
		if (is_little_endian<T>())
		{
			return in;
		}
		else
		{
			return reverse_bytes(in);
		}
	}

	template<typename T>
	inline void write_binary(std::ostream& o, const T& t)
	{
		T c = to_little_endian(t);
		o.write(reinterpret_cast<const char*>(&c), sizeof(c));
	}

	template<typename T>
	inline T read_binary(std::istream& i)
	{
		T t;
		i.read(reinterpret_cast<char*>(&t), sizeof(T));
		return to_little_endian(t);
	}

	template<typename T>
	inline void read_binary(std::istream& i, T& t)
	{
		i.read(reinterpret_cast<char*>(&t), sizeof(T));
		t = to_little_endian(t);
	}

	/**
	 * The gather operation reorders elements in the `values` array based on the indices
	 * in the `indices` array. The indices array should have length `length` and
	 * should contain integers in the range [0, `length`).
	 * After the operation, the `values` array will contain elements such that
	 * ```values[i] = previous_values[indices[i]]```
	 * The `begin` argument specifies an offset into only the `values` array.
	 * Do not add `begin` to each element in `indices`.
	 */
	template<typename T>
	void gather(std::vector<T>& values, const std::vector<size_t>& indices, size_t begin, size_t length)
	{
		std::vector<T> copy(values.begin() + begin, values.begin() + begin + length);
		for (size_t i = begin; i < begin + length; i++)
		{
			values[i] = copy[indices[i - begin]];
		}
	}

	bool load_planet_data(HostPlanetPhaseSpace& pl, const Configuration& config, std::istream& plin);
	bool load_data(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, const Configuration& config);

	size_t stable_partition_alive_indices(const std::vector<uint16_t>& flags, size_t begin, size_t length, std::unique_ptr<std::vector<size_t>>* indices);

	void save_data(const HostPlanetSnapshot& pl, const HostParticlePhaseSpace& pa, const Configuration& config, const std::string& outfile);
	void save_data_swift(const HostPlanetSnapshot& pl, const HostParticlePhaseSpace& pa, std::ostream& plout, std::ostream& icsout);

	void read_configuration(std::istream& in, Configuration* out);
	void write_configuration(std::ostream& in, const Configuration& config);

	void save_binary_track(std::ostream& trackout, const HostPlanetSnapshot& pl, const HostParticleSnapshot& pa, double time, bool to_elements, bool barycentric_elements);

	struct TrackReader
	{
		enum class State
		{
			Start,
			PlanetsBegin,
			PlanetsEnd,
			ParticlesBegin,
			ParticlesEnd,
			Finish
		};

		std::istream& input;

		inline TrackReader(std::istream& _input) : input(_input), state(State::Start) { }

		State state;

		double time;
		HostPlanetSnapshot planets;
		HostParticleSnapshot particles;

		size_t n_planets;
		size_t n_particles;

		void read_time();

		void begin_planets();
		void read_planets(const std::vector<uint32_t>* filter);
		void end_planets();

		void begin_particles();
		void read_particles();
		bool read_particle(uint32_t id);
		void end_particles();

		void check_state(const State& expected);

		static size_t bsearch_track(std::istream& f, size_t npa, uint32_t partnum, size_t stride);
	};

	struct TrackReaderOptions
	{
		bool take_all_particles;
		std::vector<uint32_t> particle_filter;
		bool take_all_planets;
		std::vector<uint32_t> planet_filter;
		double max_time;
		bool silent;

		TrackReaderOptions() : take_all_particles(false), take_all_planets(true), max_time(std::numeric_limits<double>::infinity()), silent(false) { }
	};

	void load_binary_track(std::istream& trackin, HostPlanetSnapshot& pl, HostParticleSnapshot& pa, double& time, bool skipplanets, bool skipparticles);

	void process_track(std::istream& input,
		TrackReaderOptions& options,
		const std::function<void(HostPlanetSnapshot&, HostParticleSnapshot&, double)>& callback);

	void read_tracks(const std::string& path,
		const TrackReaderOptions& options,
		const std::function<void(HostPlanetSnapshot&, HostParticleSnapshot&, double)>& callback);

	const size_t TRACK_PARTICLE_STRIDE = 28;
	const size_t TRACK_PLANET_STRIDE = 28;
}
}
