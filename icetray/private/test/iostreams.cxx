/**
   copyright (C) 2008 Troy D. Straszheim

   $Id: iostreams.cxx 168750 2016-06-01 18:26:57Z cweaver $

   @author troy d. straszheim <troy@resophonic.com>

   This file is part of IceTray.

   IceTray is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   IceTray is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with IceTray; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <I3Test.h>
#include <icetray/IcetrayFwd.h>
#include <icetray/I3Logging.h>
#include <icetray/I3Frame.h>
#include <icetray/I3Int.h>
#include <icetray/serialization.h>

#include <icetray/open.h>
#include <string>
#include <fstream>

#include <archive/binary_iarchive.hpp>
#include <archive/binary_oarchive.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
namespace io = boost::iostreams;

using namespace std;

TEST_GROUP(iostreams);

#define THROWS(CMD) try { CMD; FAIL("that should have thrown"); }	\
  catch (const std::exception &e) { /* good. */ }

TEST(zero)
{
  io::filtering_ostream fos;
  const char * file = "sink";

  io::file_sink fs(file);
  fos.push(fs);
  for (unsigned i =0; i< 100000; i++)
    fos << "HI" << i << "\n";
  fos.flush();

  remove(file);
}  

TEST(one)
{
  io::filtering_ostream fos;
  const char * file = "sink.gz";

  io::file_sink fs(file);
  fos.push(io::gzip_compressor(9));
  fos.push(fs);
  
  for (unsigned i =0; i< 100000; i++)
    fos << "HI" << i << "\n";
  fos.flush();

  remove(file);
}  

TEST(two)
{
  io::filtering_ostream fos;
  const char * file = "sinky.gz";

  I3::dataio::open(fos, file);
  
  for (unsigned i =0; i< 100000; i++)
    fos << "HI" << i << "\n";
  fos.flush();

  remove(file);
}  

TEST(three)
{
  std::string filename("framey.i3");

  I3Frame f;
  io::filtering_ostream fos;
  I3::dataio::open(fos, filename);
  f.save(fos);
  fos.flush();
  fos.reset();

  I3FramePtr rv(new I3Frame);
  io::filtering_istream fis;
  I3::dataio::open(fis, filename);
  bool r = rv->load(fis);
  ENSURE(r, "egh, failed to load frame from file");

  remove(filename.c_str());
}

TEST(i3int_only)
{
  I3IntPtr ip(new I3Int(7777));
  I3FrameObjectPtr fop(ip);
  const char * file = "i3int";

  std::ofstream ofs(file, ios::binary);
  {
    icecube::archive::portable_binary_oarchive boa(ofs);
    boa << fop;
  }
  ofs.close();
  fop.reset();
  ip.reset();

  std::ifstream ifs(file, ios::binary);
  {
    icecube::archive::portable_binary_iarchive bia(ifs);
    bia >> fop;
  }
  ip = boost::dynamic_pointer_cast<I3Int>(fop);

  ENSURE((bool)ip, "egh, not ip");
  ENSURE_EQUAL(ip->value, 7777, "egh, wrong value");

  remove(file);
}

TEST(vectorchar_only)
{
  const char * file = "vectorchar";
  vector<char> vc;
  for (char letter='A'; letter <= 'Z'; letter++)
    vc.push_back(letter);

  std::ofstream ofs(file, ios::binary);
  {
    icecube::archive::portable_binary_oarchive boa(ofs);
    boa << vc;
  }
  ofs.close();

  vector<char> vc2;
  
  std::ifstream ifs(file, ios::binary);
  {
    icecube::archive::portable_binary_iarchive bia(ifs);
    bia >> vc2;
  }
  for (int i=0; i<26; i++)
    ENSURE(vc2[i] == 'A' + i);

  remove(file);
}

TEST(six)
{
  I3Frame f;
  const char * file = "framez.i3";

  f.Put("int", I3IntPtr(new I3Int(7474)));

  std::ofstream ofs(file, ios::binary);
  f.save(ofs);
  ofs.close();

  I3Frame f2;
  std::ifstream ifs(file, ios::binary);

  f2.load(ifs);
  
  I3IntConstPtr ip = f2.Get<I3IntConstPtr>("int");

  ENSURE_EQUAL(ip->value, 7474);
  remove(file);
}

