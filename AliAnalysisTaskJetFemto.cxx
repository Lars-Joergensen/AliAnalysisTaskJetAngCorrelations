#include "AliAnalysisTaskJetFemto.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTask.h"
#include "AliPIDResponse.h"
#include "TLorentzVector.h"
#include "AliEventCuts.h"
#include "TDatabasePDG.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "THnSparse.h"
#include "AliVVZERO.h"
#include "TObjArray.h"
#include "AliAODv0.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TF1.h"
#include <vector>
#include <iostream>
#include <algorithm>


using namespace std;
ClassImp(AliAnalysisTaskJetFemto)

//__________________________________________________________________________________________________________________________________
AliAnalysisTaskJetFemto::AliAnalysisTaskJetFemto():
AliAnalysisTaskSE(),
fAODEvent(nullptr),
fPIDResponse(nullptr),
fOutputList(nullptr),
fQAList(nullptr),
fRunData(kFALSE),
hNumberOfEvents(nullptr),
hNumberOfTracks(nullptr),
hParticleIDs(nullptr),
hLeadingIDs(nullptr),
hProtonYield(nullptr),
hTPCnsigma(nullptr),
hTOFnsigma(nullptr),
hDCAxy(nullptr),
hFullPt(nullptr),
hJetPt(nullptr),
// fLeadingIDs(nullptr),
fMaximumPt(1.0),
fJetRadius(0.3)
{}
//__________________________________________________________________________________________________________________________________
AliAnalysisTaskJetFemto::AliAnalysisTaskJetFemto(const char *name):
AliAnalysisTaskSE(name),
fAODEvent(nullptr),
fPIDResponse(nullptr),
fOutputList(nullptr),
fQAList(nullptr),
fRunData(kFALSE),
hNumberOfEvents(nullptr),
hNumberOfTracks(nullptr),
hParticleIDs(nullptr),
hLeadingIDs(nullptr),
hProtonYield(nullptr),
hTPCnsigma(nullptr),
hTOFnsigma(nullptr),
hDCAxy(nullptr),
hFullPt(nullptr),
hJetPt(nullptr),
// fLeadingIDs(nullptr),
fMaximumPt(1.0),
fJetRadius(0.3)
{
    DefineInput (0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
}
//__________________________________________________________________________________________________________________________________
AliAnalysisTaskJetFemto::~AliAnalysisTaskJetFemto()  {

    fOutputList->Clear();
    delete fAODEvent;
    delete fPIDResponse;
    delete fOutputList;
    delete fQAList;
    delete hNumberOfEvents;
    delete hNumberOfTracks;
    delete hParticleIDs;
    delete hLeadingIDs;
    delete hProtonYield;
    delete hTPCnsigma;
    delete hTOFnsigma;
    delete hDCAxy;
    delete hFullPt;
    // delete fLeadingIDs;
    delete hJetPt;
}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFemto::UserCreateOutputObjects()  {

    //Create Output List
    fOutputList = new TList();
    fQAList     = new TList();
    fOutputList -> SetOwner();
    fQAList     -> SetOwner();

    //Event Counter and Centrality Distribution
    hNumberOfEvents = new TH1F ("hNumberOfEvents","hNumberOfEvents",20,0,20);
    hNumberOfEvents -> Sumw2();
    fOutputList     -> Add(hNumberOfEvents);

    hNumberOfTracks = new TH1F ("hNumberOfTracks","hNumberOfTracks",20,0,20);
    hNumberOfTracks -> Sumw2();
    fOutputList     -> Add(hNumberOfTracks);

    hParticleIDs    = new TH2F ("hParticleIDs","hParticleIDs",66001,-33000,33000, 100,0,10);
    hParticleIDs    -> Sumw2();
    fOutputList     -> Add(hParticleIDs);

    hLeadingIDs     = new TH2F ("hLeadingIDs","hLeadingIDs",66001,-33000,33000, 100,0,10);
    hLeadingIDs     -> Sumw2();
    fOutputList     -> Add(hLeadingIDs);

    hProtonYield    = new TH1F ("hProtonYield", "hProtonYield", 200,-1,1);
    hProtonYield    -> Sumw2();
    fOutputList     -> Add(hProtonYield);

    hFullPt         = new TH1D ("hFullPt", "hFullPt", 100, 0, 10);
    hFullPt         -> Sumw2();
    fOutputList     -> Add(hFullPt);

    hJetPt          = new TH1D ("hJetPt", "hJetPt", 100, 0, 10);
    hJetPt          -> Sumw2();
    fOutputList     -> Add(hJetPt);

    // add binning arrays maybe idk what to do here either
    hTPCnsigma      = new TH1F ("hTPCnsigma", "hTPCnsigma", 50,0,50);
    hTPCnsigma      -> Sumw2();
    fOutputList     -> Add(hTPCnsigma);

    hTOFnsigma      = new TH1F ("hTOFnsigma", "hTOFnsigma", 50,0,50);
    hTOFnsigma      -> Sumw2();
    fOutputList     -> Add(hTOFnsigma);

    hDCAxy          = new TH1F ("hDCAyxy", "hDCAxy", 50,0,50);
    hDCAxy          -> Sumw2();
    fOutputList     -> Add(hDCAxy);

    // fLeadingIDs     -> SetString("");
    // fOutputList     -> Add(fLeadingIDs);

    //Post Data
    PostData(1, fOutputList);
    PostData(2, fQAList);
}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFemto::UserExec(Option_t *)  {

    //Get Input Event
    if (!GetEvent()) return;

    if (fRunData) RunData ();
}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFemto::RunData()  {
    //Select tracks and put into vector. Criteria?
    AliAODEvent *fAODEvent = dynamic_cast<AliAODEvent*>(InputEvent());
    if (!fAODEvent) Fatal("AliAnalysisTaskJetFemto::RunData", "No AOD event found, check the event handler.");

    vector<Int_t> particle_ID;
    Double_t pTmax(0);
    Int_t leading_ID;
    AliAODTrack leading_particle;
    vector<AliAODTrack> protons;

    for (Int_t i=0; i<fAODEvent->GetNumberOfTracks(); i++) {
        AliAODTrack *track = static_cast<AliAODTrack*>(fAODEvent->GetTrack(i));
        if(!track) {
            continue;
            hNumberOfTracks->Fill(18.5);
        }
        hNumberOfTracks -> Fill(0.5);

        hFullPt         -> Fill(track->Pt());
        // if(!track->Charge()<=0) continue;
        hNumberOfTracks -> Fill(1.5);

        // Double_t nSigmaTPC_prot = fPIDResponse -> NumberOfSigmasTPC(track, AliPID::kProton);
        // Double_t nSigmaTOF_prot = fPIDResponse -> NumberOfSigmasTOF(track, AliPID::kProton);
        // Double_t nSigmaITS_prot = fPIDResponse -> NumberOfSigmasITS(track, AliPID::kProton);
        // // Double_t nSigmaTPC_pion = fPIDResponse -> NumberOfSigmasTPC(track, AliPID::kPion);
        // // Double_t nSigmaTOF_pion = fPIDResponse -> NumberOfSigmasTOF(track, AliPID::kPion);
        // // Double_t nSigmaITS_pion = fPIDResponse -> NumberOfSigmasITS(track, AliPID::kPion);

        // hTPCnsigma -> Fill(nSigmaTPC_prot);
        // hTOFnsigma -> Fill(nSigmaTOF_prot);
        if(TMath::Abs(track->Eta())>0.8) continue;
        hNumberOfTracks -> Fill(2.5);

        if(track->Pt()>pTmax) {
            pTmax = track->Pt();
            leading_particle = *track;
            leading_ID = i;
            hNumberOfTracks -> Fill(3.5);
            // fLeadingIDs->SetString(fLeadingIDs->GetString() + Form("%i \n", leading_ID));
        }

        if(!IsProtonCandidate(track)) continue;
        hNumberOfTracks -> Fill(4.5);
        Double_t mass = AliPID::ParticleMass(AliPID::kProton);
        hProtonYield -> Fill(GetRapidity(track, mass));

        particle_ID.push_back(i);
        protons.push_back(*track);
        // hParticleIDs->Fill(i, track->Pt());
    }

    if ((Int_t)particle_ID.size()<2) return;
    if (pTmax<fMaximumPt) return;

    // vector<AliAODTrack*> jetParticipants;
    vector<Int_t> jet_particle_ID;
    vector<AliAODTrack> jet;
    jet_particle_ID.push_back(leading_ID);
    hLeadingIDs->Fill(leading_ID, pTmax);

    Int_t exit(0);
    Int_t nPartAssociated(0);
    // Double_t pt_leading(0);

    //4-Momentum of Leading Particle
    TLorentzVector P_jet (0.,0.,0.,0.);
    P_jet.SetPxPyPzE(leading_particle.Px(),leading_particle.Py(),leading_particle.Pz(),leading_particle.E());
    Double_t pt_leading = leading_particle.Pt();

    //Jet Finder
    do {
        //Initialization
        Double_t distance_jet_min(1e+08);
        Double_t distance_bkg_min(1e+08);
        Int_t label_jet_particle(0);
        AliAODTrack jet_candidate;
        Int_t i_jet_particle(0);

        for (Int_t i=0 ; i<(Int_t)particle_ID.size() ; i++)  { //loop over all eligible particles (protons right now)
            hNumberOfTracks->Fill(5.5);
            //Skip Leading Particle & Elements already associated to the Jet
            if (particle_ID[i]==leading_ID) continue;
            if (particle_ID[i]==-1) continue;
            hNumberOfTracks->Fill(6.5);

            //Get Particle 4-Momentum
            TLorentzVector P_particle(0,0,0,0);
            AliAODTrack particle = (AliAODTrack) protons.at(i);
            P_particle.SetPxPyPzE(particle.Px(),particle.Py(),particle.Pz(),particle.E());

            //Variables
            Double_t one_over_pt2_part = 1.0/(P_particle.Pt()*P_particle.Pt());
            Double_t one_over_pt2_lead = 1.0/(P_jet.Pt()*P_jet.Pt());
            Double_t deltaY   = P_particle.Rapidity()-P_jet.Rapidity();
            Double_t deltaPhi = TVector2::Phi_0_2pi(P_particle.Phi()-P_jet.Phi());
            Double_t min      = Minimum (one_over_pt2_part,one_over_pt2_lead);
            Double_t Delta2   = deltaY*deltaY + deltaPhi*deltaPhi;

            //Distances
            Double_t distance_jet = min*Delta2/(fJetRadius*fJetRadius); // 1/pT^2 * (∆Y^2 + ∆Phi^2)/r_Jet^2
            Double_t distance_bkg = one_over_pt2_part;

            //Find Minimum Distance Jet
            if (distance_jet<distance_jet_min)  { // is current particle's dist to jet smaller than min dist
                distance_jet_min=distance_jet;// update distance between jet particles -> gets smaller
                label_jet_particle=particle_ID[i]; // remember ID
                jet_candidate = particle; // find just 1 candidate with entire for loop, the closest one
                i_jet_particle = i; // so we can remember i outside the loop
                hNumberOfTracks->Fill(7.5);
            }

            //Find Minimum Distance Bkg
            if (distance_bkg<distance_bkg_min)  { // is current particle's dist to bkg smaller than min dist
                distance_bkg_min=distance_bkg;// update min distance between jet and background particles -> gets smaller
                hNumberOfTracks->Fill(8.5);
            }
        } //for (Int_t i=0 ; i<(Int_t)particle_ID.size() ; i++)

        // Looped over relevant particles (protons at the moment) and received the next jet candidate----------------------------------------

        if (distance_jet_min<=distance_bkg_min)  {// if distance track-jet is smaller than track-background, add to jet

            //Add Particle to Jet
            jet_particle_ID.push_back(label_jet_particle);
            jet.push_back(jet_candidate);

            //Update 4-Momentum of Leading Particle
            TLorentzVector P_i(0,0,0,0);
            Double_t px_i = jet_candidate.Px();
            Double_t py_i = jet_candidate.Py();
            Double_t pz_i = jet_candidate.Pz();
            Double_t E_i  = jet_candidate.E();
            P_i.SetPxPyPzE(px_i,py_i,pz_i,E_i);
            P_jet = P_jet + P_i; //update jet momentum

            //Remove Element
            particle_ID[i_jet_particle] = -1;
            nPartAssociated++;
        } //if (distance_jet_min<=distance_bkg_min)

        if (nPartAssociated>=((Int_t)particle_ID.size()-1)) exit=1; // if counter is bigger than size of particle ID vector, exit loop - we've gone through all the particles
        if (distance_jet_min>distance_bkg_min) exit=2; // if updated jet dist is bigger than updated bkg dist, exit loop, no new jet candidate found

    } while (exit==0); // in every pass through this while loop, loop over all remaining protons, for each proton calculate dist to jet and to bkg, check if min distances need to be updated
    // after looping over all particles, min distances and a new jet candidate are found -> proceed in while loop by checking if candidate is actually part of jet
    // if it is, add to jet vector and remove from proton vector, proceed with remaining particles after resetting distances
    // continue finding candidates until run enough times or until particles are closer to background than to jet

    for (Int_t i=0; i<(Int_t)jet.size();i++) {
        hJetPt->Fill(jet.at(i).Pt());
    } //for (Int_t i=0; i<(Int_t)jet_particle_ID.size();i++)


    //Check
    if (fAODEvent->GetNumberOfTracks()==0) hNumberOfEvents->Fill(18.5);

    //Post Output Data
    PostData(1, fOutputList);
}
//__________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskJetFemto::GetEvent ()  {
    //Get Input Event
    fAODEvent = dynamic_cast <AliAODEvent*>(InputEvent());
    if (!fAODEvent) return kFALSE;
    hNumberOfEvents -> Fill(0.5);

    //Load PID Response
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();

    return kTRUE;
}
//__________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskJetFemto::IsProtonCandidate (AliAODTrack *track)  {

    //Initialization
    Bool_t isProton=(kFALSE);

    Double_t DCAxy = GetDCAtoPrimaryVertex (track,0);

    // if (!PassedTrackSelection(track,0)) return isProton;
    if (!IsHighPurityProton(track))     return isProton;
    if (TMath::Abs(DCAxy)>0.1)          return isProton;
    // if (track->P()>1.0)                 return isProton;

    isProton=kTRUE;
    return isProton;
}
//__________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskJetFemto::IsHighPurityProton (AliAODTrack *track)  {

    //Variables
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kProton);
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kProton);
    Double_t nsigmaITS = fPIDResponse -> NumberOfSigmasITS (track,AliPID::kProton);
    // Double_t nsigmaITS_recalib = GetRecalibratedITSnsigma (nsigmaITS,track->Eta(),track->P());
    Double_t pt = track->Pt();

    //Selection
    if (pt<0.7 /*&& TMath::Abs(nsigmaITS_recalib)<3.0*/ && TMath::Abs(nsigmaTPC)<3.0) return kTRUE;
    if (pt>0.7 && TMath::Abs(nsigmaTPC)<2.0 && TMath::Abs(nsigmaTOF)<2.0)         return kTRUE;

    return kFALSE;
}
//__________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskJetFemto::GetRapidity (AliAODTrack *track, Double_t mass)  {

    //Initialization
    Double_t y(-999);

    //Rapidity Calculation
    Double_t p  = track->P();
    Double_t pz = track->Pz();
    Double_t E = TMath::Sqrt(mass*mass + p*p);
    if (E == TMath::Abs(pz)) return -999;
    y = 0.5*TMath::Log((E+pz)/(E-pz));

    return y;
}
//_______________________________________________________________________________________________________________________________________
Double_t  AliAnalysisTaskJetFemto::Minimum (Double_t x1, Double_t x2)  {

    Double_t x_min(x1);
    if (x1<x2) x_min = x1;
    if (x1>x2) x_min = x2;

    return x_min;
}
//__________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskJetFemto::GetDCAtoPrimaryVertex (AliAODTrack *track, Int_t index)  {

    Double_t dca[2];
    Double_t covMatrix[3];
    if (!track->PropagateToDCA (fAODEvent->GetPrimaryVertex(),fAODEvent->GetMagneticField(),10000,dca,covMatrix)) return -999;
    hDCAxy->Fill(dca[0]);

    return dca[index];
}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFemto::Terminate(Option_t *)  {

    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
//__________________________________________________________________________________________________________________________________
