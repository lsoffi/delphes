/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/** \class PileUpMerger
 *
 *  Merges particles from pile-up sample into event
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "modules/PileUpMerger.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesTF2.h"
#include "classes/DelphesPileUpReader.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>



using namespace std;

//------------------------------------------------------------------------------

PileUpMerger::PileUpMerger() :
  fFunction(0), fReader(0), fItInputArray(0)
{
  fFunction = new DelphesTF2;
}


//------------------------------------------------------------------------------

PileUpMerger::~PileUpMerger()
{
  delete fFunction;
}

//------------------------------------------------------------------------------

void PileUpMerger::Init()
{
  const char *fileName;

  fTimeResolution = GetDouble("TimeResolution", 30E-12);
  fTimeWindow = GetDouble("TimeWindow", 10);
  fLowEtaRange = GetDouble("LowEtaRange", 1.5);
  fHighEtaRange = GetDouble("HighEtaRange", 1.7);
  fEfficiencyTiming = GetDouble("EfficiencyTiming", 1.0);

  fPileUpDistribution = GetInt("PileUpDistribution", 0);

  fMeanPileUp  = GetDouble("MeanPileUp", 10);
  
  fRadius = GetDouble("Radius", 1.0);
  fRadius2 = fRadius*fRadius;
  fHalfLength = GetDouble("HalfLength", 3.0);
  fBz = GetDouble("Bz", 0.0);
  if(fRadius < 1.0E-2)
  {
    cout << "ERROR: magnetic field radius is too low\n";
    return;
  }
  if(fHalfLength < 1.0E-2)
  {
    cout << "ERROR: magnetic field length is too low\n";
    return;
  }

  fRadiusMax = GetDouble("RadiusMax", fRadius);
  fHalfLengthMax = GetDouble("HalfLengthMax", fHalfLength);


  fZVertexSpread = GetDouble("ZVertexSpread", 0.15);
  fTVertexSpread = GetDouble("TVertexSpread", 1.5E-09);

  fInputBeamSpotX = GetDouble("InputBeamSpotX", 0.0);
  fInputBeamSpotY = GetDouble("InputBeamSpotY", 0.0);
  fOutputBeamSpotX = GetDouble("OutputBeamSpotX", 0.0);
  fOutputBeamSpotY = GetDouble("OutputBeamSpotY", 0.0);

  // read vertex smearing formula

  fFunction->Compile(GetString("VertexDistributionFormula", "0.0"));
  fFunction->SetRange(-fZVertexSpread, -fTVertexSpread, fZVertexSpread, fTVertexSpread);

  fileName = GetString("PileUpFile", "MinBias.pileup");
  fReader = new DelphesPileUpReader(fileName);

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Delphes/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();
  // import input array
  fInputArraySecond = ImportArray(GetString("InputArraySecond", "Delphes/stableParticles"));
  fItInputArraySecond = fInputArray->MakeIterator();

  // create output arrays
  fParticleOutputArray = ExportArray(GetString("ParticleOutputArray", "stableParticles"));
  fVertexOutputArray = ExportArray(GetString("VertexOutputArray", "vertices"));
}

//------------------------------------------------------------------------------

void PileUpMerger::Finish()
{
  if(fReader) delete fReader;
}

//------------------------------------------------------------------------------

void PileUpMerger::Process()
{
  TDatabasePDG *pdg = TDatabasePDG::Instance();
  TParticlePDG *pdgParticle;
  Int_t pid, nch, nvtx = -1;
  Float_t x, y, z, t, vx, vy, vt, vt_ref;
  Float_t px, py, pz, e, pt;
  Double_t dz, dphi, dt, sumpt2, dz0, dt0;
  Int_t numberOfEvents, event, numberOfParticles;
  Long64_t allEntries, entry;
  Candidate *candidate, *candidateSecond, *vertex;
  DelphesFactory *factory;

  const Double_t c_light = 2.99792458E8;

  fItInputArray->Reset();

  // --- Deal with primary vertex first  ------

  fFunction->GetRandom2(dz, dt);

  dz0 = -1.0e6;
  dt0 = -1.0e6;

  dt *= c_light*1.0E3; // necessary in order to make t in mm/c
  dz *= 1.0E3; // necessary in order to make z in mm

  vx = 0.0;
  vy = 0.0;
  vt = 0.0;

  numberOfParticles = fInputArray->GetEntriesFast();
  nch = 0;
  sumpt2 = 0.0;

  factory = GetFactory();
  vertex = factory->NewCandidate();

  std::cout << "I come here often " << numberOfParticles << std::endl;

  //while((candidate = static_cast<Candidate*>(fItInputArray->Next())) && (candidateSecond = static_cast<Candidate*>(fItInputArraySecond->Next())))
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    std::cout << "but not here " << numberOfParticles << std::endl;
    vx += candidate->Position.X();
    vy += candidate->Position.Y();
    z = candidate->Position.Z();
    t = candidate->Position.T();
    pt = candidate->Momentum.Pt();
  

    //===========================Original===============================
    //if (dz0 < -999999.0)
    //  dz0 = z;
    //if (dt0 < -999999.0)
    //  dt0 = t;
    //// cancel any possible offset in position and time the input file
    //candidate->Position.SetZ(z - dz0 + dz);
    //candidate->Position.SetT(t - dt0 + dt);
    //=================================================================
    
    //=======================New Added by Pablo========================
    Double_t tf_smeared = gRandom->Gaus(0, fTimeResolution);
    dt = dt + tf_smeared * 1.0E3*c_light;
    candidate->Position.SetZ(z + dz);
    candidate->Position.SetT(t + dt);
    candidate->ErrorT = 200E-12 * 1.0E3 * c_light;
    //Propagate the particle to the timing detector geometry
    //Candidate *theProp = static_cast<Candidate*>(candidate->Clone());
    //Candidate *prop = Propagate(theProp);
    //if(prop != NULL) {
    //    Double_t eta = fabs(prop->Momentum.Eta());
    //    //If the particle is not in the eta gap of the timing detector or it was detected given the efficiency we update the error 
    //    if(!(eta > fLowEtaRange && eta < fHighEtaRange) && (gRandom->Uniform(1.0) < fEfficiencyTiming)) { 
    //        candidate->ErrorT = fTimeResolution * 1.0E3 * c_light;
    //    }
    //}
    vt += (t + dt);
    candidate->IsPU = 0;
    fParticleOutputArray->Add(candidate);
    std::cout << "Particle put in array" << std::endl;
    if(TMath::Abs(candidate->Charge) >  1.0E-9)
    {
      //Added by Pablo
      nch++;
      sumpt2 += pt*pt;
      vertex->AddCandidate(candidate);
    }
  }


  if(nch > 0)
  {
    vx /= nch;
    vy /= nch;
    vt /= nch;
  }

  nvtx++;
  vertex->Position.SetXYZT(vx, vy, dz, vt);
  vertex->ClusterIndex = nvtx;
  vertex->ClusterNDF = nch;
  vertex->SumPT2 = sumpt2;
  vertex->GenSumPT2 = sumpt2;
  fVertexOutputArray->Add(vertex);

  //Added by Pablo: contains the time of the primary vertex.
  vt_ref = vt;
  
  // --- Then with pile-up vertices  ------

  switch(fPileUpDistribution)
  {
    case 0:
      numberOfEvents = gRandom->Poisson(fMeanPileUp);
      break;
    case 1:
      numberOfEvents = gRandom->Integer(2*fMeanPileUp + 1);
      break;
    case 2:
      numberOfEvents = fMeanPileUp;
      break;
    default:
      numberOfEvents = gRandom->Poisson(fMeanPileUp);
      break;
  }

  allEntries = fReader->GetEntries();


  for(event = 0; event < numberOfEvents; ++event)
  {
    do
    {
      entry = TMath::Nint(gRandom->Rndm()*allEntries);
    }
    while(entry >= allEntries);

    fReader->ReadEntry(entry);

   // --- Pile-up vertex smearing

    fFunction->GetRandom2(dz, dt);

    dt *= c_light*1.0E3; // necessary in order to make t in mm/c
    dz *= 1.0E3; // necessary in order to make z in mm

    dphi = gRandom->Uniform(-TMath::Pi(), TMath::Pi());

    vx = 0.0;
    vy = 0.0;

    numberOfParticles = 0;
    sumpt2 = 0.0;
    nch = 0;

    //factory = GetFactory();
    vertex = factory->NewCandidate();
    while(fReader->ReadParticle(pid, x, y, z, t, px, py, pz, e))
    {
      candidate = factory->NewCandidate();

      candidate->PID = pid;

      candidate->Status = 1;

      pdgParticle = pdg->GetParticle(pid);
      candidate->Charge = pdgParticle ? Int_t(pdgParticle->Charge()/3.0) : -999;
      candidate->Mass = pdgParticle ? pdgParticle->Mass() : -999.9;

      //=======================New Added by Pablo========================
      Double_t tf_smeared = gRandom->Gaus(0, fTimeResolution);
      dt = dt + tf_smeared * 1.0E3*c_light;
      vt += (t + dt);
 
      candidate->IsPU = 1;

      candidate->Momentum.SetPxPyPzE(px, py, pz, e);
      candidate->Momentum.RotateZ(dphi);
      
      Candidate *theProp = static_cast<Candidate*>(candidate->Clone());
      Candidate *extParticle = Propagate(theProp);
      pt = candidate->Momentum.Pt();

      x -= fInputBeamSpotX;
      y -= fInputBeamSpotY;
      candidate->Position.SetXYZT(x, y, z + dz, t + dt);
      candidate->ErrorT = 200E-12 * 1.0E3 * c_light;
      //Propagate the particle to the timing detector geometry
      //Candidate *prop = Propagate(candidate);
      //Bool_t timeMeasurement = false;
      //if(prop != NULL) { 
      //    Double_t eta = fabs(prop->Momentum.Eta());
      //    //If the particle is not in the eta gap of the timing detector or it was detected given the efficiency we update the error 
      //    if(!(eta > fLowEtaRange && eta < fHighEtaRange) && (gRandom->Uniform(1.0) < fEfficiencyTiming)) { 
      //      candidate->ErrorT = fTimeResolution * 1.0E3 * c_light;
      //      timeMeasurement = true;
      //    }
      //}
       
      candidate->Position.RotateZ(dphi);
      candidate->Position += TLorentzVector(fOutputBeamSpotX, fOutputBeamSpotY, 0.0, 0.0);

      vx += candidate->Position.X();
      vy += candidate->Position.Y();
      //If it's a charged particle, with a valid timeMeasurement, and outside a time window w.r.t. vertex
      //if(timeMeasurement && TMath::Abs(candidate->Charge) >  1.0E-9 && TMath::Abs((t + dt) - vt_ref) > fTimeWindow * 1.0E3 * c_light) continue; 

      ++numberOfParticles;
      if(TMath::Abs(candidate->Charge) >  1.0E-9)
      {
        nch++;
        sumpt2 += pt*pt;
        vertex->AddCandidate(candidate);
      }

      fParticleOutputArray->Add(candidate);
    }
    //Added by Pablo: don't want to consider vertices where there are no charged particles. 
    if(nch == 0) continue;
    if(nch > 0)
    {
      vx /= nch;
      vy /= nch;
      vt /= nch;
    }

    nvtx++;

    vertex->Position.SetXYZT(vx, vy, dz, dt);

    vertex->ClusterIndex = nvtx;
    vertex->ClusterNDF = nch;
    vertex->SumPT2 = sumpt2;
    vertex->GenSumPT2 = sumpt2;

    vertex->IsPU = 1;

    fVertexOutputArray->Add(vertex);

  }
}



Candidate *PileUpMerger::Propagate(Candidate *_particle)
{
  Candidate *candidate, *mother, *particle;
  TLorentzVector particlePosition, particleMomentum, beamSpotPosition;
  Double_t px, py, pz, pt, pt2, e, q;
  Double_t x, y, z, t, r, phi;
  Double_t x_c, y_c, r_c, phi_c, phi_0;
  Double_t x_t, y_t, z_t, r_t;
  Double_t t1, t2, t3, t4, t5, t6;
  Double_t t_z, t_r, t_ra, t_rb;
  Double_t tmp, discr, discr2;
  Double_t delta, gammam, omega, asinrho;
  Double_t rcu, rc2, xd, yd, zd;
  Double_t l, d0, dz, p, ctgTheta, phip, etap, alpha;
  Double_t bsx, bsy, bsz;

  const Double_t c_light = 2.99792458E8;

  candidate = _particle;
  particle = _particle;
  particlePosition = particle->Position;
  particleMomentum = particle->Momentum;


  x = particlePosition.X()*1.0E-3;
  y = particlePosition.Y()*1.0E-3;
  z = particlePosition.Z()*1.0E-3;

  bsx = beamSpotPosition.X()*1.0E-3;
  bsy = beamSpotPosition.Y()*1.0E-3;
  bsz = beamSpotPosition.Z()*1.0E-3;

  q = particle->Charge;

  // check that particle position is inside the cylinder
  if(TMath::Hypot(x, y) > fRadiusMax || TMath::Abs(z) > fHalfLengthMax)
  {
    return NULL;
  }

  px = particleMomentum.Px();
  py = particleMomentum.Py();
  pz = particleMomentum.Pz();
  pt = particleMomentum.Pt();
  pt2 = particleMomentum.Perp2();
  e = particleMomentum.E();

  if(pt2 < 1.0E-9)
  {
    return NULL;
  }

  if(TMath::Hypot(x, y) > fRadius || TMath::Abs(z) > fHalfLength)
  {
    mother = candidate;
    candidate = static_cast<Candidate*>(candidate->Clone());
    candidate->InitialPosition = particlePosition;
    candidate->Position = particlePosition;
    candidate->L = 0.0;

    candidate->Momentum = particleMomentum;
    candidate->AddCandidate(mother);
    return candidate;
  }
  else if(TMath::Abs(q) < 1.0E-9 || TMath::Abs(fBz) < 1.0E-9)
  {
    // solve pt2*t^2 + 2*(px*x + py*y)*t - (fRadius2 - x*x - y*y) = 0
    tmp = px*y - py*x;
    discr2 = pt2*fRadius2 - tmp*tmp;

    if(discr2 < 0.0)
    {
      // no solutions
      return NULL;
    }

    tmp = px*x + py*y;
    discr = TMath::Sqrt(discr2);
    t1 = (-tmp + discr)/pt2;
    t2 = (-tmp - discr)/pt2;
    t = (t1 < 0.0) ? t2 : t1;

    z_t = z + pz*t;
    if(TMath::Abs(z_t) > fHalfLength)
	    {
      t3 = (+fHalfLength - z) / pz;
      t4 = (-fHalfLength - z) / pz;
      t = (t3 < 0.0) ? t4 : t3;
    }

    x_t = x + px*t;
    y_t = y + py*t;
    z_t = z + pz*t;

    l = TMath::Sqrt( (x_t - x)*(x_t - x) + (y_t - y)*(y_t - y) + (z_t - z)*(z_t - z));

    mother = candidate;
    candidate = static_cast<Candidate*>(candidate->Clone());

    candidate->InitialPosition = particlePosition;
    candidate->Position.SetXYZT(x_t*1.0E3, y_t*1.0E3, z_t*1.0E3, particlePosition.T() + t*e*1.0E3);
    candidate->L = l*1.0E3;

    candidate->Momentum = particleMomentum;
    candidate->AddCandidate(mother);

    return candidate;
  }
  else
  {

    // 1.  initial transverse momentum p_{T0}: Part->pt
    //     initial transverse momentum direction phi_0 = -atan(p_X0/p_Y0)
    //     relativistic gamma: gamma = E/mc^2; gammam = gamma * m
    //     gyration frequency omega = q/(gamma m) fBz
    //     helix radius r = p_{T0} / (omega gamma m)

    gammam = e*1.0E9 / (c_light*c_light);      // gammam in [eV/c^2]
    omega = q * fBz / (gammam);                // omega is here in [89875518/s]
    r = pt / (q * fBz) * 1.0E9/c_light;        // in [m]

    phi_0 = TMath::ATan2(py, px); // [rad] in [-pi, pi]

    // 2. helix axis coordinates
    x_c = x + r*TMath::Sin(phi_0);
    y_c = y - r*TMath::Cos(phi_0);
    r_c = TMath::Hypot(x_c, y_c);
    phi_c = TMath::ATan2(y_c, x_c);
    phi = phi_c;
    if(x_c < 0.0) phi += TMath::Pi();

    rcu = TMath::Abs(r);
    rc2 = r_c*r_c;

    // calculate coordinates of closest approach to track circle in transverse plane xd, yd, zd
    xd = x_c*x_c*x_c - x_c*rcu*r_c + x_c*y_c*y_c;
    xd = (rc2 > 0.0) ? xd / rc2 : -999;
    yd = y_c*(-rcu*r_c + rc2);
    yd = (rc2 > 0.0) ? yd / rc2 : -999;
    zd = z + (TMath::Sqrt(xd*xd + yd*yd) - TMath::Sqrt(x*x + y*y))*pz/pt;

    // use perigee momentum rather than original particle
    // momentum, since the orignal particle momentum isn't known

    px = TMath::Sign(1.0, r) * pt * (-y_c / r_c);
    py = TMath::Sign(1.0, r) * pt * (x_c / r_c);
    etap = particleMomentum.Eta();
    phip = TMath::ATan2(py, px);

    particleMomentum.SetPtEtaPhiE(pt, etap, phip, particleMomentum.E());

    // calculate additional track parameters (correct for beamspot position)

    d0        = ((x - bsx) * py - (y - bsy) * px) / pt;
    dz        = z - ((x - bsx) * px + (y - bsy) * py) / pt * (pz / pt);
    p         = particleMomentum.P();
    ctgTheta  = 1.0 / TMath::Tan (particleMomentum.Theta());


    // 3. time evaluation t = TMath::Min(t_r, t_z)
    //    t_r : time to exit from the sides
    //    t_z : time to exit from the front or the back
    t_r = 0.0; // in [ns]
    int sign_pz = (pz > 0.0) ? 1 : -1;
    if(pz == 0.0) t_z = 1.0E99;
    else t_z = gammam / (pz*1.0E9/c_light) * (-z + fHalfLength*sign_pz);

    if(r_c + TMath::Abs(r)  < fRadius)
    {
      // helix does not cross the cylinder sides
      t = t_z;
    }
    else
    {
      asinrho = TMath::ASin((fRadius*fRadius - r_c*r_c - r*r) / (2*TMath::Abs(r)*r_c));
      delta = phi_0 - phi;
      if(delta <-TMath::Pi()) delta += 2*TMath::Pi();
      if(delta > TMath::Pi()) delta -= 2*TMath::Pi();
      t1 = (delta + asinrho) / omega;
      t2 = (delta + TMath::Pi() - asinrho) / omega;
      t3 = (delta + TMath::Pi() + asinrho) / omega;
      t4 = (delta - asinrho) / omega;
      t5 = (delta - TMath::Pi() - asinrho) / omega;
      t6 = (delta - TMath::Pi() + asinrho) / omega;

      if(t1 < 0.0) t1 = 1.0E99;
      if(t2 < 0.0) t2 = 1.0E99;
      if(t3 < 0.0) t3 = 1.0E99;
      if(t4 < 0.0) t4 = 1.0E99;
      if(t5 < 0.0) t5 = 1.0E99;
      if(t6 < 0.0) t6 = 1.0E99;

      t_ra = TMath::Min(t1, TMath::Min(t2, t3));
      t_rb = TMath::Min(t4, TMath::Min(t5, t6));
      t_r = TMath::Min(t_ra, t_rb);
      t = TMath::Min(t_r, t_z);
    }

    // 4. position in terms of x(t), y(t), z(t)
    x_t = x_c + r * TMath::Sin(omega * t - phi_0);
    y_t = y_c + r * TMath::Cos(omega * t - phi_0);
    z_t = z + pz*1.0E9 / c_light / gammam * t;
    r_t = TMath::Hypot(x_t, y_t);


    // compute path length for an helix

    alpha = pz*1.0E9 / c_light / gammam;
    l = t * TMath::Sqrt(alpha*alpha + r*r*omega*omega);

    if(r_t > 0.0)
    {

      // store these variables before cloning
      if(particle == candidate)
      {
        particle->D0 = d0*1.0E3;
        particle->DZ = dz*1.0E3;
        particle->P = p;
        particle->PT = pt;
        particle->CtgTheta = ctgTheta;
        particle->Phi = phip;
      }

      mother = candidate;
      candidate = static_cast<Candidate*>(candidate->Clone());

      candidate->InitialPosition = particlePosition;
      candidate->Position.SetXYZT(x_t*1.0E3, y_t*1.0E3, z_t*1.0E3, particlePosition.T() + t*c_light*1.0E3);

      candidate->Momentum = particleMomentum;

      candidate->L = l*1.0E3;

      candidate->Xd = xd*1.0E3;
      candidate->Yd = yd*1.0E3;
      candidate->Zd = zd*1.0E3;

      candidate->AddCandidate(mother);

      return candidate;
    }
  }
}


