#include "AliAnalysisTaskJetFemto.h"
#include "AliInputEventHandler.h" // maybe
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h" // maybe
#include "AliAnalysisTask.h"
#include "TLorentzVector.h"
#include "AliPIDResponse.h"
#include "AliEventCuts.h" // maybe
#include "TDatabasePDG.h"
//#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "TObjArray.h" // maybe
#include "TVector2.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TF1.h"
#include <vector>
#include <iostream>
#include <algorithm>


using namespace std;
ClassImp(AliAnalysisTaskJetFemto)

AliAnalysisTaskJetFemto::AliAnalysisTaskJetFemto(): AliAnalysisTaskSE(),
    fESDevent(nullptr),
    fOutputList(nullptr),
    fQAList(nullptr),
    fPIDResponse(nullptr),
    hNumberOfEvents(nullptr),
    hTestHistPt(nullptr)
    {}
//_______________________________________________________________________________________________________________________________________
AliAnalysisTaskJetFemto::AliAnalysisTaskJetFemto(const char *name):
    AliAnalysisTaskSE(name),
    fESDevent(nullptr),
    fOutputList(nullptr),
    fQAList(nullptr),
    fPIDResponse(nullptr),
    hNumberOfEvents(nullptr),
    hTestHistPt(nullptr)
    {
        DefineInput (0, TChain::Class());
        DefineOutput(1, TList::Class());
        DefineOutput(2, TList::Class());
    }
//_______________________________________________________________________________________________________________________________________
AliAnalysisTaskJetFemto::~AliAnalysisTaskJetFemto()  {

    fOutputList->Clear();
    delete fESDevent;
    delete fOutputList;
    delete fQAList;
    delete fPIDResponse;
    delete hTestHistPt;
    delete hNumberOfEvents;
}
//_______________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFemto::UserCreateOutputObjects()  {

    //Create Output List
    fOutputList = new TList();
    fQAList     = new TList();
    fOutputList -> SetOwner();
    fQAList     -> SetOwner();
    hTestHistPt = new TH1D("hTestHist", "hTestHist", 100, 0, 100);
    hTestHistPt -> Sumw2();
    hNumberOfEvents = new TH1D("hNumberOfEvents", "hNumberOfEvents", 1, 0, 1);
    hNumberOfEvents -> Sumw2();
    fOutputList -> Add(hTestHistPt);
    fOutputList -> Add(hNumberOfEvents);

    //Post Data
    PostData(1, fOutputList);
    PostData(2, fQAList);
}
//_______________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFemto::UserExec(Option_t *)  {
    fESDevent = dynamic_cast<AliESDEvent*>(InputEvent());
    if(!fESDevent) ::Fatal("AliAnalysisTaskJetFemto::UserExec", "No ESD event found, check the event handler.");

    int iTracks{fESDevent->GetNumberOfTracks()};
    for(int i{0}; i < iTracks; i++) {
        AliESDtrack* track = static_cast<AliESDtrack*>(fESDevent->GetTrack(i)); // get i-th track and put it into the variable "track"
        if(!track) continue;

        //Double_t nsigmaITS = fPIDResponse -> NumberOfSigmasITS (track,AliPID::kProton);
        //Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kProton);
        //Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kProton);

        if(track->Pt()>2.0) continue;
        if(track->Charge()<0) continue;
        //if(nsigmaTPC()<3) continue;
        hTestHistPt->Fill(track->Pt());
        hNumberOfEvents->Fill(0.5);
    }
    //Jet Finder
    //NEW OBJECTS: distance_jet_min, distance_bkg_min, label_jet_particle, i_jet_particle, particle_ID, leading_ID, P_Particle,
    /*do {

        //Initialization
        Double_t distance_jet_min(1e+08);
        Double_t distance_bkg_min(1e+08);
        Int_t label_jet_particle(0);
        Int_t i_jet_particle(0);

        for (Int_t i=0 ; i<(Int_t)particle_ID.size() ; i++)  {

            //Skip Leading Particle & Elements already associated to the Jet
            if (particle_ID[i]==leading_ID) continue;
            if (particle_ID[i]==-1)         continue;

            //Get Particle 4-Momentum
            TLorentzVector P_particle(0,0,0,0);
            AliMCParticle *particle = (AliMCParticle*) fMCEvent->GetTrack(particle_ID[i]);
            P_particle.SetPxPyPzE(particle->Px(),particle->Py(),particle->Pz(),particle->E());

            //Variables
            Double_t one_over_pt2_part = 1.0/(P_particle.Pt()*P_particle.Pt());
            Double_t one_over_pt2_lead = 1.0/(P_leading.Pt()*P_leading.Pt());
            Double_t deltaY   = P_particle.Rapidity()-P_leading.Rapidity();
            Double_t deltaPhi = TVector2::Phi_0_2pi(P_particle.Phi()-P_leading.Phi());
            Double_t min      = Minimum (one_over_pt2_part,one_over_pt2_lead);
            Double_t Delta2   = deltaY*deltaY + deltaPhi*deltaPhi;

            //Distances
            Double_t distance_jet = min*Delta2/(fJetRadius*fJetRadius);
            Double_t distance_bkg = one_over_pt2_part;

            //Find Minimum Distance Jet
            if (distance_jet<distance_jet_min)  {
                distance_jet_min=distance_jet;
                label_jet_particle=particle_ID[i];
                i_jet_particle = i;
            }

            //Find Minimum Distance Bkg
            if (distance_bkg<distance_bkg_min)  {
                distance_bkg_min=distance_bkg;
            }
        }

        if (distance_jet_min<=distance_bkg_min)  {

            //Add Particle to Jet
            jet_particle_ID.push_back(label_jet_particle);

            //Update 4-Momentum of Leading Particle
            TLorentzVector P_i(0,0,0,0);
            AliMCParticle *particle = (AliMCParticle*) fMCEvent->GetTrack(label_jet_particle);
            Double_t px_i = particle->Px();
            Double_t py_i = particle->Py();
            Double_t pz_i = particle->Pz();
            Double_t E_i  = particle->E();
            P_i.SetPxPyPzE(px_i,py_i,pz_i,E_i);
            P_leading = P_leading + P_i;

            //Remove Element
            particle_ID[i_jet_particle] = -1;
            nPartAssociated++;
        }

        if (nPartAssociated>=((Int_t)particle_ID.size()-1)) exit=1;
        if (distance_jet_min>distance_bkg_min) exit=2;

    } while (exit==0);*/
}
//_______________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFemto::Terminate(Option_t *)  {
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}