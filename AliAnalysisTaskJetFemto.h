#ifndef AliAnalysisTaskJetFemto_cxx
#define AliAnalysisTaskJetFemto_cxx

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTask.h"
#include "TLorentzVector.h"
#include "AliESDEvent.h"
#include "TVector3.h"
#include "TList.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TF1.h"

class AliAnalysisTaskJetFemto : public AliAnalysisTaskSE {
public:
   AliAnalysisTaskJetFemto();
   AliAnalysisTaskJetFemto(const char *name);
   virtual ~AliAnalysisTaskJetFemto();

   // General Functions
   virtual void UserCreateOutputObjects();
   virtual void UserExec(Option_t *option);
   virtual void Terminate(Option_t *);

private: // exclamation marks mean transient members, none means persistent --> ROOT streamer
   AliESDEvent *fESDevent;   //!
   TList       *fOutputList; //!
   TList       *fQAList;     //!
   Double_t     fMaximumPt;  //
   Double_t     fJetRadius;  //
   AliPIDResponse  *fPIDResponse;//!

   //TF1 *fProtonWeights; //

   // Event Counter
   TH1D *hNumberOfEvents; //!

   // p_{T} Spectra
   TH1D *hTestHistPt;                  //!

   AliAnalysisTaskJetFemto(const AliAnalysisTaskJetFemto &);
   AliAnalysisTaskJetFemto &operator=(const AliAnalysisTaskJetFemto &);

   ClassDef(AliAnalysisTaskJetFemto, 1);
};

#endif