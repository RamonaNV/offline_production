
#include "corsika-reader/CorsikaClient.h"

#include <boost/asio.hpp>
#include <boost/process.hpp>

#include <iostream>

using boost::asio::ip::tcp;
namespace fs = boost::filesystem;
namespace process = boost::process;

class CorsikaClientImpl {
public:
	CorsikaClientImpl() : acceptor_(io_, tcp::endpoint(tcp::v4(), 0)), socket_(io_)
	{}
	boost::asio::io_service io_;
	boost::asio::ip::tcp::acceptor acceptor_;
	boost::asio::ip::tcp::socket socket_;
	boost::process::child corsika_;

};

namespace {

void configure(CorsikaClientImpl &impl, const CorsikaClient::config_type &config, boost::process::opstream &steering)
{
  steering
    << "RUNNR " << 17 << "\n"      // run number
    << "EVTNR " << 1 << "\n"       // number of first shower event
    << "NSHOW " << INT_MAX << "\n" // number of showers to generate
    << "PRMPAR 14\n"               // primary type (not used)
    << "ESLOPE -1\n"               // slope of primary energy spectrum (not used)
    << "ERANGE 1.E6  1.E6\n"       // energy range of primary particle (not used)
    << "THETAP  0.  0.\n"          // range of zenith angle (degree, not used)
    << "PHIP -180.  180.\n"        //range of azimuth angle (degree, not used)
    ;
  for (auto &stream : config.rng) {
    steering << "SEED " << stream.seed
             << " " << stream.ncalls%1000000000
             << " " << stream.ncalls/1000000000 << "\n";
  }
  steering
    << "OBSLEV  283199\n"           // observation level (in cm)
    << "ELMFLG  F   T    \n"        // em. interaction flags (NKG,EGS)
    //ARRANG
    //FIXHEI
    //FIXCHI
    << "MAGNET  16.4  -53.4\n"      // magnetic field south pole
    << "HADFLG  0  1  0  1  0  2\n" // flags hadr.interact.&fragmentation
    //SIBYLL
    //SIBSIG
    ;
  /*
  {
    steering << "ECUTS ";
    for (int i=0; i < 4; i++){
      steering << config.elcuts[i] << " ";      
    }
    steering << "\n";
    }*/

  steering << "ECUTS   273.0000 273.0000 0.0030 0.0030\n";
  
  steering
    << "*MUADDI  T        \n"       // additional info for muons
    << "*NUADDI  T       \n"        // additional info for neutrinos
    << "MUMULT  T        \n"        // muon multiple scattering angle
    << "LONGI F 1.E10 F F\n"        // longit.distr. & step size & fit & out
    << "MAXPRT  0        \n"        // max. number of printed events
    << "ECTMAP  1.E21    \n"        // cut on gamma factor for printout
    << "STEPFC  1.0      \n"        // mult. scattering step length fact.
    << "DEBUG   F  6  F  1000000\n" // debug flag and log.unit for out
    //PIPE
    //DETCFG
    << "ATMOD " << config.atmosphere << "\n" // atmospheric model
    << "DIRECT  ./       \n"        // output directory
    
    //options not in simprod card
    << "RADNKG  200.E2   \n"        // outer radius for NKG lat.dens.distr.
    << "DATBAS  F        \n"        // write .dbase file
    << "PAROUT  F F      \n"        // suppress DAT file
    << "USER    you      \n"        // user 
    << "DYNSTACK 1024    \n"
    << "DYNSTACK_P bias_target mu bias_factor " << 1 << "\n"
    << "REMOTE_CONTROL 127.0.0.1 " << impl.acceptor_.local_endpoint().port() << "\n"
    << "EXIT" << std::endl
    ;
  // send an EOF to start argument processing
  steering.pipe().close();
}


enum message_id : uint32_t {
	UNKNOWN = 0x00,
	EVENT_HEADER = 0x03,
	EVENT_END = 0x04,
	RUN_HEADER = 0x05,
	RUN_END = 0x06,
	PARTICLE = 0x07,
	PRIMARY = 0x08,
};

// RemoteControl packet:
// - Packet length (4 bytes)
// - Packet ID (4 bytes)
// - Payload (Packet length - 8 bytes)
struct packet_header {
	packet_header() : length(0), type(UNKNOWN) {};
	uint32_t length;
	message_id type;
};

packet_header
peek_message(boost::asio::ip::tcp::socket &socket)
{
	packet_header header;
	socket.read_some(boost::asio::buffer(&header, 8));
	assert(header.length >= 8);
	header.length -= 8;
	return header;
}

template <typename T>
std::vector<T>
read(boost::asio::ip::tcp::socket &socket, uint32_t bytes)
{
	assert(bytes % sizeof(T) == 0);
	std::vector<T> buffer(bytes/sizeof(T));
	socket.read_some(boost::asio::buffer(buffer.data(), bytes));
	return buffer;
}

void
seek(boost::asio::ip::tcp::socket &socket, uint32_t bytes)
{
	std::vector<char> buffer(bytes);
	socket.read_some(boost::asio::buffer(buffer.data(), bytes));
}

};

CorsikaClient::config_type::config_type() :
    atmosphere(11), elcuts({{320.,320.,1e21,1e21}}),
    rng({{{1,0},{2,0}}})
{}

CorsikaClient::CorsikaClient(const std::string &corsika_executable, config_type config) :
    impl_(new CorsikaClientImpl),
    more_particles_available_(false)
{
	process::opstream steering;
	impl_->corsika_ = process::child(corsika_executable,
	    process::std_out > process::null, process::std_err > stderr,
	    // process::std_out > stdout, process::std_err > stderr,
	    process::std_in < steering,
	    process::start_dir(fs::path(corsika_executable).parent_path()));
	configure(*impl_, config, steering);
	
	impl_->acceptor_.accept(impl_->socket_);
}

CorsikaClient::~CorsikaClient() {}

std::vector<double>
CorsikaClient::StartShower(uint32_t particle_id, double energy, double theta,
                           double phi, I3ShowerBias::BiasParticleType bias_target,
                           double bias_factor, std::array<double,4> elcuts)
{
	{
		std::array<double,6> initial = {{double(particle_id),energy,theta,phi,double(bias_target),bias_factor}};
		std::array<uint32_t,2> header = {{ 8+sizeof(initial)+sizeof(elcuts), PRIMARY }};
		std::array<boost::asio::const_buffer,3> bufs = {{
			boost::asio::buffer(header),
			boost::asio::buffer(initial),
			boost::asio::buffer(elcuts)
		}};
		impl_->socket_.send(bufs);
	}
	
	event_header_.clear();
	do {
		packet_header next_packet = peek_message(impl_->socket_);
		if (next_packet.type != EVENT_HEADER) {
			seek(impl_->socket_, next_packet.length);
		} else {
			event_header_ = read<float>(impl_->socket_, next_packet.length);
		}
	} while (event_header_.empty());
	
	std::vector<double> primary;
	{
		packet_header next_packet = peek_message(impl_->socket_);
		assert(next_packet.type == PRIMARY);
		primary = read<double>(impl_->socket_, next_packet.length);
	}
	
	event_trailer_.clear();
	more_particles_available_ = true;
	
	return primary;
}

std::vector<double>
CorsikaClient::NextParticle()
{
	std::vector<double> particle;

	while (particle.empty() && more_particles_available_) {
		packet_header next_packet = peek_message(impl_->socket_);
		if (next_packet.type == PARTICLE) {
			particle = read<double>(impl_->socket_, next_packet.length);
		} else if (next_packet.type == EVENT_END) {
			event_trailer_ = read<float>(impl_->socket_, next_packet.length);
			more_particles_available_ = false;
		} else {
			seek(impl_->socket_, next_packet.length);
		}
	} 
	
	return particle;
}

const std::vector<float>&
CorsikaClient::GetEventHeader() const
{
	return event_header_;
}

const std::vector<float>&
CorsikaClient::GetEventEnd() const
{
	return event_trailer_;
}

