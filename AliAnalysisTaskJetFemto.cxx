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
fJetRadius(0.5)
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
fJetRadius(0.5)
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

    hParticleIDs    = new TH1I ("hParticleIDs","hParticleIDs",66001,-33000,33000);
    hParticleIDs    -> Sumw2();
    fOutputList     -> Add(hParticleIDs);

    hLeadingIDs    = new TH1I ("hLeadingIDs","hLeadingIDs",66001,-33000,33000);
    hLeadingIDs    -> Sumw2();
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

    // leadingParticles = new vector<AliAODTrack*>;

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

    for (Int_t i=0; i<fAODEvent->GetNumberOfTracks(); i++) {
        AliAODTrack *track = static_cast<AliAODTrack*>(fAODEvent->GetTrack(i));
        if(!track) {
            continue;
            hNumberOfTracks->Fill(17.5);
        }
        hNumberOfTracks -> Fill(0.5);

        hFullPt         -> Fill(track->Pt());
        if(track->Charge()<=0) continue;
        hNumberOfTracks -> Fill(1.5);

        Double_t nSigmaTPC_prot = fPIDResponse -> NumberOfSigmasTPC(track, AliPID::kProton);
        Double_t nSigmaTOF_prot = fPIDResponse -> NumberOfSigmasTOF(track, AliPID::kProton);
        Double_t nSigmaITS_prot = fPIDResponse -> NumberOfSigmasITS(track, AliPID::kProton);
        Double_t nSigmaTPC_pion = fPIDResponse -> NumberOfSigmasTPC(track, AliPID::kPion);
        Double_t nSigmaTOF_pion = fPIDResponse -> NumberOfSigmasTOF(track, AliPID::kPion);
        Double_t nSigmaITS_pion = fPIDResponse -> NumberOfSigmasITS(track, AliPID::kPion);

        hTPCnsigma -> Fill(nSigmaTPC_prot);
        hTOFnsigma -> Fill(nSigmaTOF_prot);
        if(TMath::Abs(track->Eta())>0.8) continue;
        hNumberOfTracks -> Fill(2.5);

        if(track->Pt()>2.0 && track->Pt()>pTmax) {
            pTmax = track->Pt();
            leading_ID = track->GetID();
            hNumberOfTracks -> Fill(3.5);
        }

        if(!IsProtonCandidate(track)) continue;
        hNumberOfTracks -> Fill(4.5);
        Double_t mass = AliPID::ParticleMass(AliPID::kProton);
        hProtonYield -> Fill(GetRapidity(track, mass));

        particle_ID.push_back(track->GetID());
        hParticleIDs->Fill(track->GetID());
    }

    // vector<AliAODTrack*> jetParticipants;
    vector<Int_t> jet_particle_ID;
    jet_particle_ID.push_back(leading_ID);
    hLeadingIDs->Fill(leading_ID);
    Int_t exit(0);
    Int_t nPartAssociated(0);

    //4-Momentum of Leading Particle
    // TLorentzVector P_leading (0,0,0,0);
    // AliAODTrack *leading_particle = (AliAODTrack*) fAODEvent->GetTrack(leading_ID); //take leading_ID and create leading particle TLorentzVec
    // if(!leading_particle) {
    //     return;
    //     hNumberOfEvents->Fill(1.5);
    // }
    // if(leading_particle) {hNumberOfEvents->Fill(2.5);}
    // P_leading.SetPxPyPzE(leading_particle->Px(),leading_particle->Py(),leading_particle->Pz(),leading_particle->E());
    // Double_t pt_leading = P_leading.Pt();

    //Jet Finder
    // do {
    //     //Initialization
    //     Double_t distance_jet_min(1e+08);
    //     Double_t distance_bkg_min(1e+08);
    //     Int_t label_jet_particle(0);
    //     Int_t i_jet_particle(0);

        /* for (Int_t i=0 ; i<(Int_t)particle_ID.size() ; i++)  { //loop over all eligible particles (protons right now)

            //Skip Leading Particle & Elements already associated to the Jet
            if (particle_ID[i]==leading_ID) continue;
            if (particle_ID[i]==-1)         continue;

            //Get Particle 4-Momentum
            TLorentzVector P_particle(0,0,0,0);
            AliAODTrack *particle = (AliAODTrack*) fAODEvent->GetTrack(particle_ID[i]); //casts the track from the MC event into the particle variable
            P_particle.SetPxPyPzE(particle->Px(),particle->Py(),particle->Pz(),particle->E());

            //Variables
            Double_t one_over_pt2_part = 1.0/(P_particle.Pt()*P_particle.Pt());
            Double_t one_over_pt2_lead = 1.0/(P_leading.Pt()*P_leading.Pt());
            Double_t deltaY   = P_particle.Rapidity()-P_leading.Rapidity();
            Double_t deltaPhi = TVector2::Phi_0_2pi(P_particle.Phi()-P_leading.Phi());
            Double_t min      = Minimum (one_over_pt2_part,one_over_pt2_lead);
            Double_t Delta2   = deltaY*deltaY + deltaPhi*deltaPhi;

            //Distances
            Double_t distance_jet = min*Delta2/(fJetRadius*fJetRadius); // 1/pT^2 * (∆Y^2 + ∆Phi^2)/r_Jet^2
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
        } */

        // exit=1;

        /* if (distance_jet_min<=distance_bkg_min)  {

            //Add Particle to Jet
            jet_particle_ID.push_back(label_jet_particle);

            //Update 4-Momentum of Leading Particle
            TLorentzVector P_i(0,0,0,0);
            AliAODTrack *particle = (AliAODTrack*) fAODEvent->GetTrack(label_jet_particle);
            Double_t px_i = particle->Px();
            Double_t py_i = particle->Py();
            Double_t pz_i = particle->Pz();
            Double_t E_i  = particle->E();
            P_i.SetPxPyPzE(px_i,py_i,pz_i,E_i);
            P_leading = P_leading + P_i;

            //Remove Element
            particle_ID[i_jet_particle] = -1;
            nPartAssociated++;
        } */

        // if (nPartAssociated>=((Int_t)particle_ID.size()-1)) exit=1;
        // if (distance_jet_min>distance_bkg_min) exit=2;

    // } while (exit==0);

    // for(Int_t i=0; i<(Int_t)jet_particle_ID.size();i++) {
    //     AliAODTrack *particle = (AliAODTrack*) fAODEvent->GetTrack(jet_particle_ID[i]);
    //     hJetPt->Fill(particle->Pt());
    // }


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
    if (track->P()>1.0)                 return isProton;

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

    return dca[index];
}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetFemto::Terminate(Option_t *)  {

    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
//__________________________________________________________________________________________________________________________________
