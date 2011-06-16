//-*-C++-*-
// $Id: TG4BLStream.cc 149 2008-07-02 06:33:39Z comet $

static const char* TG4BLStream_cc =
"@(#) $Id: TG4BLStream.cc 149 2008-07-02 06:33:39Z comet $";

#include <iostream>
#include <iomanip>

#include "G4EventManager.hh"

#include "TG4BLStream.hh"

TG4BLStream::TG4BLStream()
  : fFileStatus(0), fFileName("primary.gun"),
    fIsValid(false),
    fPosition(G4ThreeVector(0,0,0)),
    fDirection(G4ThreeVector(0,0,0)),
    fMomentum(0), fProperTime(0), fWeight(0),
    fPDGEncoding(0), fEventID(0), fTrackID(0), fParentID(0)
{
}



TG4BLStream::TG4BLStream(G4String filename)
  : fFileStatus(0), fFileName(filename),
    fIsValid(false),
    fPosition(G4ThreeVector(0,0,0)),
    fDirection(G4ThreeVector(0,0,0)),
    fMomentum(0), fProperTime(0), fWeight(0),
    fPDGEncoding(0), fEventID(0), fTrackID(0), fParentID(0)
{
  this->Open(filename);
}



TG4BLStream::~TG4BLStream()
{
}



void TG4BLStream::Next()
{
#define BUFSIZE 500
  const G4int bufsize = BUFSIZE;
  char buffer[BUFSIZE], buff2[BUFSIZE];


  if (!fFileIn.is_open()) {
    G4cout << "TG4BLStream: File was not opened." << G4endl;
    fIsValid = false;
    return;
  }

//   G4cout << fFileIn.bad() << G4endl;

  G4int ieof = 0;
  while (true) {

    if (fFileIn.eof()) {
      if (ieof>0) {
	G4cout << "TG4BLStream: No valid data lines at all." << G4endl;
	fIsValid = false;
	break;
      }
      fFileIn.clear(); // Reset EOF status before seekg.
      fFileIn.seekg(0, std::ios::beg);
      ++ieof;
      ++fFileStatus;
      if (G4EventManager::GetEventManager()->GetVerboseLevel() > 0) {
	G4cout << "G4BLStream: an input file was rewound." << G4endl;
      }
    }

    char c, ctmp;
    c = fFileIn.get();  // get the first charactor.
    if ( c != '#' ) {
      // data line was found.
      fFileIn.putback(c);
    } else {
      // comment line.
      fFileIn.getline(buffer,bufsize);  // discard the rest of line.
      G4cout << "#" << buffer << G4endl;
      continue;
    }

    // Get a data line.
    G4double pos[3], mom[3], time, weight;
    G4int pdgid, iev, trackid, parentid;
    fFileIn.getline(buffer,bufsize);
    G4int nitem = sscanf(buffer,"%lf %lf %lf %lf %lf %lf %lf %500c",
			 &pos[0], &pos[1], &pos[2], &mom[0], &mom[1], &mom[2],
			 &time, buff2);
    G4int nitem2 = sscanf(buff2,"%d %d %d %d %lf",
			  &pdgid, &iev, &trackid, &parentid, &weight);
    if (nitem+nitem2<13) {
      // less number of items.  Probably an empty line.
      continue;
    }

    fPosition = G4ThreeVector(pos[0], pos[1], pos[2])*mm;
    G4ThreeVector momentum = G4ThreeVector(mom[0], mom[1], mom[2])*MeV;
    fDirection = momentum.unit();
    fMomentum = momentum.mag();
    fProperTime = time*ns;
    fPDGEncoding = pdgid;
    fEventID = iev;
    fTrackID = trackid;
    fParentID = parentid;
    fWeight = weight;

    fIsValid = true;
    break;
  }

//   G4cout << "Position = " << fPosition << G4endl;
//   G4cout << "Direction = " << fDirection << G4endl;
//   G4cout << "momentum = " << fMomentum << G4endl;
//   G4cout << "time = " << fProperTime << G4endl;
//   G4cout << "pdgid = " << fPDGEncoding << G4endl;
//   G4cout << "weight = " << fWeight << G4endl;
//   G4cout << "eventID, trackID, parentID = "
// 	 << fEventID << ", " << fTrackID << ", " << fParentID << G4endl;

}



void TG4BLStream::Open(G4String val)
{
  fFileName = val;

  fFileIn.open(fFileName, std::ios::in);
//   G4cout << "TG4BLStream: Open the input file : " << fFileName
// 	 << G4endl;
  if ( !fFileIn.is_open() ) {
    G4cout << "TG4BLStream: unable to open the input file : " << fFileName
	   << G4endl;
    exit(1);
  }
  fFileStatus = 1;
}



void TG4BLStream::Close()
{
  fFileIn.close();
  fFileStatus = 0;
}

