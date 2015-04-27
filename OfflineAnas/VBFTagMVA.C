#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom2.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include "THStack.h"
//#include "setTDRStyle.C"
#include "FormatMC_VBFTag.h"
#include "Cuts.h"

#define nsam 3 

#define inputSkimedVBF200 49923 
#define inputSkimedGGH200 299979 
#define inputSkimedDYJM50 30461028  
#define inputSkimedDYJM50Fraction 3.23023e5

//pb
#define X_VBF200 0.8685 
#define X_GGH200 7.127  

// use higgs-> 2l2q number  
//#define BrHZZ 2.55e-01
#define Br2l2q 3.60e-02

// X_Br= X*Br2l2q
#define X_BrVBF200 0.0313   
#define X_BrGGH200 0.2566 
#define X_BrDYJM50 3503.71  
     

float BrXsecDYJ =3503.71;
float BrXsecGGH200= 0.2566;
float BrXsecVBF200= 0.0313;

float NDYJ =30461028.0;
float NGGH200= 299979.0;
float NVBF200= 49923.0;

float weight[nsam]={0.0};

/*
struct sample_t {
char filename[128];
char tag[128];
}; 


struct sample_t  sample[nsam] = {
//{"/home/cdxfe/CDXFE_hzz2l2q/VBF2012/forTMVA/DYJetsToLL_M50_VBF_VBFTree_v1_1_1_Ijm.root","DYJ1"},
//{"/home/cdxfe/CDXFE_hzz2l2q/VBF2012/forTMVA/DYJetsToLL_M50_VBF_VBFTree_v1_1_10.root","DYJ1"},
//{"/home/cdxfe/CDXFE_hzz2l2q/VBF2012/forTMVA/VBF200_VBFTree_v1_1_4.root","VBF"},
//{"/home/cdxfe/CDXFE_hzz2l2q/VBF2012/forTMVA/GGH200_VBFTree_v1_1_2.root","GGH"}

//{"/home/cdxfe/CDXFE_hzz2l2q/VBF2012/forTMVA_v2/VBFTree_DYJ.root","DYJ"},
//{"/home/cdxfe/CDXFE_hzz2l2q/VBF2012/forTMVA_v2/VBFTree_VBF.root","VBF"},
//{"/home/cdxfe/CDXFE_hzz2l2q/VBF2012/forTMVA_v2/VBFTree_GGH.root","GGH"}

{"/home/cdxfe/CDXFE_hzz2l2q/VBF2012/forTMVA_v3/VBFTree_VBF.root","VBF"},
{"/home/cdxfe/CDXFE_hzz2l2q/VBF2012/forTMVA_v3/VBFTree_DYJ.root","DYJ"},
{"/home/cdxfe/CDXFE_hzz2l2q/VBF2012/forTMVA_v3/VBFTree_GGH.root","GGH"},

//{"/data4/yunju/TreeMC2/TreeMCPho80_v6_40to80.root","TreeMCPho80_v6_40to80.root"},
//{"/data4/yunju/TreeMC2/TreeMCPho80_v6_80to120.root","TreeMCPho80_v6_80to120.root"},
//{"/data4/yunju/TreeMC2/TreeMCPho80_v6_120to160.root","TreeMCPho80_v6_120to160.root"},
//{"/data4/yunju/TreeMC2/TreeMCPho80_v6_160to200.root","TreeMCPho80_v6_160to200.root"},
//{"/data4/yunju/TreeMC2/TreeMCPho80_v6_200to240.root","TreeMCPho80_v6_200to240.root"},
};

*/


double deltaPhi(double Phi1, double Phi2){
   double result = Phi1-Phi2;
   while (result > TMath::Pi()) result -= 2*TMath::Pi();
   while (result <= -TMath::Pi()) result += 2*TMath::Pi();
   return result;                                  }//deltaPhi


void VBFTagMVA()
{
   
   float weightVBF=0.0;float weightGGH=0.0;float weightDYJ=0.0;
   weightVBF= X_BrVBF200/inputSkimedVBF200;
   weightGGH= X_BrGGH200/inputSkimedGGH200; 
   weightDYJ= X_BrDYJM50/inputSkimedDYJM50; 

   struct sample_t {
   char filename[128];
   char tag[128];
   float weight;
   }; 


struct sample_t  sample[nsam] = {
//{"/home/cdxfe/CDXFE_hzz2l2q/VBF2012/forTMVA/DYJetsToLL_M50_VBF_VBFTree_v1_1_1_Ijm.root","DYJ1"},
//{"/home/cdxfe/CDXFE_hzz2l2q/VBF2012/forTMVA/DYJetsToLL_M50_VBF_VBFTree_v1_1_10.root","DYJ1"},
//{"/home/cdxfe/CDXFE_hzz2l2q/VBF2012/forTMVA/VBF200_VBFTree_v1_1_4.root","VBF"},
//{"/home/cdxfe/CDXFE_hzz2l2q/VBF2012/forTMVA/GGH200_VBFTree_v1_1_2.root","GGH"}

//{"/home/cdxfe/CDXFE_hzz2l2q/VBF2012/forTMVA_v3/DYJetsToLL_M50_VBF.root","DYJ",weightDYJ},
//{"/home/cdxfe/CDXFE_hzz2l2q/VBF2012/forTMVA_v3/VBF200_VBF.root","VBF",weightVBF},
//{"/home/cdxfe/CDXFE_hzz2l2q/VBF2012/forTMVA_v3/GGH200_VBF.root","GGH",weightGGH}

//{"/data4/yunju/VBF2012/v2Feb6/DYJetsToLL_M50_VBF_all.root","DYJ",weightDYJ},
//{"/data4/yunju/VBF2012/v2Feb6/VBF200_VBF.root","VBF",weightVBF},
//{"/data4/yunju/VBF2012/v2Feb6/GGH200_VBF.root","GGH",weightGGH}
{"/data4/yunju/VBF2012/v2Feb6/DYJetsToLL_M50_VBF_all.root","DYJ",1.150e-04},
{"/data4/yunju/VBF2012/v2Feb6/VBF200_VBF.root","VBF",6.26e-4},
{"/data4/yunju/VBF2012/v2Feb6/GGH200_VBF.root","GGH",8.5539e-7}


//{"/data4/yunju/TreeMC2/TreeMCPho80_v6_40to80.root","TreeMCPho80_v6_40to80.root"},
//{"/data4/yunju/TreeMC2/TreeMCPho80_v6_80to120.root","TreeMCPho80_v6_80to120.root"},
//{"/data4/yunju/TreeMC2/TreeMCPho80_v6_120to160.root","TreeMCPho80_v6_120to160.root"},
//{"/data4/yunju/TreeMC2/TreeMCPho80_v6_160to200.root","TreeMCPho80_v6_160to200.root"},
//{"/data4/yunju/TreeMC2/TreeMCPho80_v6_200to240.root","TreeMCPho80_v6_200to240.root"},
};






   weight[0] = (BrXsecDYJ/NDYJ);
   weight[1] = BrXsecVBF200/NVBF200;
   weight[2] = BrXsecGGH200/NGGH200;  

    
    
 
     gROOT->SetStyle("Plain");
   gStyle->SetOptStat(0);
     //setTDRStyle();
     //save tree for TMVA
     TTree *tree[nsam];
     TTree *Skimtree[nsam];

     FILE *ffISO[nsam];
 

     TH1F *h_sieieEB =new TH1F("sieieEB","sieieEB_",42,0,0.03);

     TH1F *h_EtaDY =new TH1F("EtaDY ","DY_",50,0,1000);
     TH1F *h_EtaVBF =new TH1F("EtaVBF ","VBF_",50,0,1000);
     TH1F *h_EtaGGH =new TH1F("EtaGGH ","GGH_",50,0,1000);

     TH1F *h_MaxEta[3]; 
     TH1F *h_MaxMjj[3]; 
     TH1F *h_mcU3[3];
     TH1F *h_mcU2[3];
     TH1F *h_mcU2eta[3];
  
     
     for(int file=0; file<nsam; file++)
     {
        h_mcU2[file] =new TH1F(Form("mcU2_%s",sample[file].tag),Form("PmcU2_%s",sample[file].tag),600,-300,300);
        h_mcU2eta[file]= new TH1F(Form("mcU2eta_%s",sample[file].tag),Form("PmcU2eta_%s",sample[file].tag),600,-300,300);
        h_MaxEta[file]= new TH1F(Form("MaxEta_%s",sample[file].tag),Form("PMaxEta_%s",sample[file].tag),50,0,8);
        h_MaxMjj[file]= new TH1F(Form("MaxMjj_%s",sample[file].tag),Form("PMaxMjj_%s",sample[file].tag),50,0,1000);

        if(file>runNMC) continue;       
   	//if(strcmp(sample[file].tag,"May10rereco30")==0) continue;	
        TFile *f1 = new TFile(sample[file].filename);
        
        ffISO[file] = fopen(Form("MVA_%s_v3.dat",sample[file].tag),"w");
        fprintf(ffISO[file],"var1/F:var2/F\n");
   //      TFile *f1 = new TFile("/data4/yunju/TreeMCPho120/TreeMCPho120_v6_0to40.root");
  //     TFile *f1 = new TFile("/home/yunju/DataTran/PhotreeOct10/TreeMCPho80_v4/TreeMCPho80_v4.root");
  //     TFile *f1 = new TFile("/home/yunju/DataTran/PhotreeOct10/TreeMCPho50_v4/TreeMCPho50_v4.root");
        //for TMVA
        float Mjj=0.0;
        float etajj=0.0;  
        float SkimMjj=0.0; 
        int nMeta=0;
        int nMeta1tag=0;
        int nMeta0tag=0;
        int nMMjj=0;
        int nMMjj1tag=0;
        int nMMjj0tag=0;
        int nTa=0;
        int totEvt=0;
        int MatchedEvt=0;
        tree[file] = new TTree(Form("MVA_%s",sample[file].tag),Form("MVA_%s",sample[file].tag));
        tree[file]->Branch("Mjj_",&Mjj,"Mjj/F");
        tree[file]->Branch("etajj_",&etajj,"etajj/F");
        Skimtree[file] = new TTree(Form("SkimMVA_%s",sample[file].tag),"test_tree2");
        Skimtree[file]->Branch("SkimMjj",&SkimMjj,"SkimMjj/F");


        TTree *root = new TTree();
        root = (TTree*)f1->FindObjectAny("tree");

        // create the varibles and  branchs so that we can  read the tree the from file 
	EvtInfoBranches Inf;
        Inf.Register(root);
        printf("%s ; weight : %f  ;%d entries to be loaded \n",sample[file].tag,sample[file].weight ,root->GetEntries());
        for(Long64_t jentry=0;jentry<root->GetEntriesFast();jentry++)
        { 
        //     if(jentry>2000) continue;
            //for root tree
            root->GetEntry(jentry);
            // trigger selection

            //cout<<" Event "<< jentry <<endl;  
       /* 
            for(int n=0;n<Inf.genJetE_->size();n++)
             {
              cout<<" genJet "<<Inf.genJetE_->at(n)<<" "<<Inf.genJetE_->at(n) <<" "<<Inf.genJetPt_->at(n)<<" "<<Inf.genJetEta_->at(n)<<" "<<Inf.genJetPhi_->at(n)<<endl;  

             }
*/


            int nGoodEle=0; 
            for(int n=0;n<Inf.patElecEt_->size();n++)
            {
        
/*
       CanEle *myEle;
               myEle->SetMemberValue
               (
                 patElecPt_,
                 patElecEta_,
                 patElecDelEtaIn_,
                 patElecDelPhiIn_,
                 patElecSigIhIh_,
                 patElecHoE_,
                            ,
                            ,  
                            ,
                            ,

  
*patElecEta_;
*patElecPhi_;
*patElecM_;
*patElecScEn_;
*patElecScEt_;
*patElecScEta_;
*patElecScPhi_;
*patElecisEcalDriven_
*patElecisTrackerDriv
*patElecSigIhIh_;
*patElecDelEtaIn_;
*patElecDelPhiIn_;
*patElecHoE_;
*patElecTrkIso_;
*patElecHcalIso_;
*patElecEcalIso_;
*patElecCharge_;
*patElecRelIsoComb_;

               );


*/
               if(fabs(Inf.patElecEta_->at(n))>2.5) continue;
               if(fabs(Inf.patElecPt_->at(n))<20) continue; 
               if(fabs(Inf.patElecEta_->at(n))<1.44 && Inf.patElecHoE_->at(n)>0.12 ) continue;
               if(fabs(Inf.patElecEta_->at(n))>1.44 && Inf.patElecHoE_->at(n)>0.1 ) continue; 
               if(fabs(Inf.patElecEta_->at(n))<1.44 && Inf.patElecDelEtaIn_->at(n)>0.007 ) continue;
               if(fabs(Inf.patElecEta_->at(n))<1.44 && Inf.patElecDelPhiIn_->at(n)>0.15 ) continue;
               if(fabs(Inf.patElecEta_->at(n))<1.44 && Inf.patElecSigIhIh_->at(n)>0.01) continue;
               if(Inf.patElecID_->at(n)==0) continue;
               nGoodEle++;
            }

            int nGoodMu=0;
            
            for(int n=0;n<Inf.patMuonEta_->size();n++)
            {
               
               if(fabs(Inf.patMuonEta_->at(n))>2.4) continue;
               if(fabs(Inf.patMuonPt_->at(n))<20) continue; 
               if(Inf.patMuonID_->at(n)==0) continue;
               nGoodMu++;
            }

            if(nGoodMu < 2 && nGoodEle < 2) continue;


      
//             cout<<"Number of jets "<< Inf.JetPt_->size() <<endl;
             float MaxMjj=-100.0;
             float Maxetajj=-100.0;  
             bool foundGoodjetpair=false;
             bool PassCentralMjj=false;
             int nGoodJets=0;
             //variable for MC check
             float JetMMjj1U3=10000;
             float JetMMjj2U3=10000;
             float JetMMjj1U2=10000;
             float JetMMjj2U2=10000;
             float JetMMjj1U1=10000;
             float JetMMjj2U1=10000;
             int JetMMjj1matched =-100;   
             int JetMMjj2matched =-100;

             float JetMeta1U3=10000;
             float JetMeta2U3=10000;
             float JetMeta1U2=10000;
             float JetMeta2U2=10000;
             float JetMeta1U1=10000;
             float JetMeta2U1=10000;
             int JetMeta1matched =-100;   
             int JetMeta2matched =-100;
                 
         // cout<<"Mjj sss  :"<< Mjj<<endl;
     
        for(int n=0;n<Inf.JetPt_->size();n++)
        {
           if(Inf.JetPt_->at(n)<30) continue;
           if(fabs(Inf.JetEta_->at(n))>4.7) continue;

            //jet id 

            if(Inf.JetNeutEmEFr_->at(n)>0.99) continue;
            if(Inf.JetNeutHadEFr_->at(n)>0.99) continue;
            if(Inf.JetCharEmEFr_->at(n)>0.99 &&fabs(Inf.JetEta_->at(n))<2.4) continue; 
            if(Inf.JetCharHadEFr_->at(n)<= 0 &&fabs(Inf.JetEta_->at(n))<2.4) continue;
            if(Inf.JetCharMulti_->at(n) <= 0 &&fabs(Inf.JetEta_->at(n))<2.4) continue;

            nGoodJets++; 
            //cout<<n <<" th jet "<< Inf.JetPt_->at(n)<<" "<<Inf.JetEta_->at(n)<<" "<<Inf.JetPhi_->at(n)<<" "<<Inf.JetEn_->at(n) <<endl;
           TLorentzVector LVJetn;     
           LVJetn.SetPtEtaPhiE(Inf.JetPt_->at(n),Inf.JetEta_->at(n),Inf.JetPhi_->at(n),Inf.JetEn_->at(n));  
           for(int i=0;i<Inf.JetPt_->size();i++)
           {
              if(Inf.JetPt_->at(i)<30) continue;
              if(fabs(Inf.JetEta_->at(i))>4.7) continue;
              if(n>=i) continue;
             

              if(Inf.JetNeutEmEFr_->at(i)>0.99) continue;
              if(Inf.JetNeutHadEFr_->at(i)>0.99) continue;
              if(Inf.JetCharEmEFr_->at(i)>0.99&&fabs(Inf.JetEta_->at(i))<2.4) continue;
              if(Inf.JetCharHadEFr_->at(i)<=0 &&fabs(Inf.JetEta_->at(i))<2.4) continue;
              if(Inf.JetCharMulti_->at(i)<=0 &&fabs(Inf.JetEta_->at(i))<2.4) continue;
              
              foundGoodjetpair=true;
              TLorentzVector LVJeti;
              LVJeti.SetPtEtaPhiE(Inf.JetPt_->at(i),Inf.JetEta_->at(i),Inf.JetPhi_->at(i),Inf.JetEn_->at(i));
              TLorentzVector LVJJ; 
              LVJJ = LVJetn+LVJeti;
              //cout<<"Mjj Mine  :"<< LVJJ.M()<<endl;
               Mjj=LVJJ.M();
               if(LVJJ.M()>75&&LVJJ.M()<105)
               {PassCentralMjj=true;}
               //cout<<"Mjjmmmmvhvhm sss  :"<< Mjj<<endl;   
               if(LVJJ.M()>MaxMjj) 
               {
           
                 MaxMjj=LVJJ.M(); 
                 JetMMjj1matched = Inf.JethasGenParton_->at(i);
                 JetMMjj2matched = Inf.JethasGenParton_->at(n);   
                 JetMMjj1U3=Inf.JetGenPartonU3ID_->at(i);
                 JetMMjj2U3=Inf.JetGenPartonU3ID_->at(n);
                 JetMMjj1U2=Inf.JetGenPartonU2ID_->at(i);
                 JetMMjj2U2=Inf.JetGenPartonU2ID_->at(n);
                 JetMMjj1U1=Inf.JetGenPartonU1ID_->at(i);
                 JetMMjj2U1=Inf.JetGenPartonU1ID_->at(n);
             
               }
              // cout<<"Mjjmmmmm sss  :"<< Mjj<<endl;
               etajj=fabs(Inf.JetEta_->at(i)-Inf.JetEta_->at(n));
              // if(file==0)h_EtaDY->Fill(Maxetajj);
              // if(file==1)h_EtaVBF->Fill(Maxetajj);
              //  if(file==2)h_EtaGGH->Fill(Maxetajj);   
               if(etajj>Maxetajj)   
               {
                 Maxetajj=  etajj;
                
                 JetMeta1matched = Inf.JethasGenParton_->at(i);
                 JetMeta2matched = Inf.JethasGenParton_->at(n);   
              
                 JetMeta1U3=Inf.JetGenPartonU3ID_->at(i);
                 
                 JetMeta2U3=Inf.JetGenPartonU3ID_->at(n);
                 JetMeta1U2=Inf.JetGenPartonU2ID_->at(i);
                 JetMeta2U2=Inf.JetGenPartonU2ID_->at(n);
                 JetMeta1U1=Inf.JetGenPartonU1ID_->at(i);
                 JetMeta2U1=Inf.JetGenPartonU1ID_->at(n);
                 
                }

               //tree[file]->Fill();
         
          }//jet i
      

      }//jet n

/*
cout<<"Mjeti: "<<JetMMjj1U1<<" "<<JetMMjj1U2<<" "<<JetMMjj1U3<<endl;
cout<<"Mjetj "<<JetMMjj2U1<<" "<<JetMMjj2U2<<" "<<JetMMjj2U3<<endl;
cout<<"Etajeti: "<<JetMeta1U1<<" "<<JetMeta1U2<<" "<<JetMeta1U3<<endl;
cout<<"Etajetj: "<<JetMeta2U1<<" "<<JetMeta2U2<<" "<<JetMeta2U3<<endl;
*/
//if(foundGoodjetpair&&JetMMjj1U1!=-99999&&JetMMjj2U1!=-99999&&JetMeta1U1!=-99999&&JetMeta2U1!=-99999) MatchedEvt++; 



//cin.get();
if(foundGoodjetpair&&JetMeta1matched&&JetMeta2matched&&JetMMjj1matched&&JetMMjj2matched){
  MatchedEvt++;
  if(JetMMjj1U3==2212&&JetMMjj2U3==2212&&JetMeta1U3==2212&&JetMeta2U3==2212) nTa++;
  if(JetMMjj1U3==2212&&JetMMjj2U3==2212)
  { 
     nMMjj++;
  }
  else if(JetMMjj1U3!=2212&&JetMMjj2U3!=2212) 
  { 
    nMMjj0tag++;
    //cout<<JetMMjj1U1<<"  "<<JetMMjj1U2<<" "<<JetMMjj1U3<<endl;
    h_mcU2[file]->Fill(JetMMjj1U2);
    h_mcU2[file]->Fill(JetMMjj2U2);
    //cin.get();
   }
   else 
   {
    nMMjj1tag++;
   }

  if(JetMeta1U3==2212&&JetMeta2U3==2212) 
  {
    nMeta++;
  }
  else if(JetMeta1U3!=2212&&JetMeta2U3!=2212)
  {
     nMeta0tag++;
     h_mcU2eta[file]->Fill(JetMeta1U2);
     h_mcU2eta[file]->Fill(JetMeta2U2);
  }
  else
  {
    nMeta1tag++;
  }

}

       if(foundGoodjetpair&&PassCentralMjj) 
       {
       fprintf(ffISO[file],"%f %f\n",MaxMjj,Maxetajj);
       h_MaxEta[file]->Fill(Maxetajj);
       h_MaxMjj[file]->Fill(MaxMjj);
       }
           // tree[file]->Fill(); 
      for(int n=0;n<Inf.higgsPt->size();n++)
      {
          //cout<<"Official :"<< Inf.zjjM->at(n)<<endl;
          SkimMjj=Inf.zjjM->at(n);
              //Skimtree[file]->Fill();  
      }
      for(int n=0;n<3;n++)
      {
         SkimMjj=1000;
 //        Skimtree[file]->Fill();
      }
 
      
      for(int j=0;j<Inf.HjetIndex->size();j++)
      { 
            //  cout<<"Jetindex :"<<Inf.HjetIndex->at(j)<<" "<< Inf.HjetPt->at(j)  <<" "<<Inf.HjetEta->at(j)<<" "<<Inf.HjetPhi->at(j)<<" "<<Inf.HjetE->at(j) <<endl;
              
      }   
      if(jentry%10000==1)
      cout<<"##### Event "<< jentry <<endl;
//         cin.get();  
        
     totEvt=jentry;
    }// end of looping over entries
    fclose(ffISO[file]);
    cout<<" Sum "<< totEvt<<" " <<"; Matched :  "<<MatchedEvt<<"; "<<nTa<<" "<<nMeta0tag<<" "<<nMeta1tag<<" "<<nMeta<<" "<<nMMjj<<endl;
    cout<<" Sum "<< totEvt<<" " <<"; Matched :  "<<(float)MatchedEvt/totEvt<<"; "<<(float)nTa/MatchedEvt<<" "<<(float)nMeta0tag/MatchedEvt<<" "<<(float)nMeta1tag/MatchedEvt<<" "<<(float)nMeta/MatchedEvt<<" Mjj : "<<(float)nMMjj0tag/MatchedEvt<<" "<<(float)nMMjj1tag/MatchedEvt<<" "<<(float)nMMjj/MatchedEvt<<endl;
 //   cin.get();
        
  }//file loop
    //gStyle->SetPalette(1);
  TFile*f = new TFile("VBFKinVar.root","recreate"); 
  for(int n=0;n<nsam;n++)
  {

   if(n>runNMC) continue;
//   tree[n]->Write();
//   tree[n]->Delete();
//   Skimtree[n]->Write();
//   Skimtree[n]->Delete();
    
    h_MaxEta[n]->Write();
    h_MaxMjj[n]->Write();

   //cout<<"weight :"<< 1e9*weight[n]<<endl;
 // h_mcU2[n]->Write(); 
 // h_mcU2eta[n]->Write();

    

  }
    //h_EtaDY->Write();
    //h_EtaVBF->Write();
    //h_EtaGGH->Write();   

  f->Close();

 TLegend *leg = new TLegend(0.65, 0.65, 0.90, 0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.03);


TH1F *h_SMaxEta[3];
TH1F *h_SMaxMjj[3];  
for(int n=0;n<nsam;n++)
{ 
  if(n>runNMC) continue;
  h_SMaxEta[n]=(TH1F*)h_MaxEta[n]->Clone();
  h_SMaxEta[n]->Scale(sample[n].weight);
  h_SMaxMjj[n]=(TH1F*)h_MaxMjj[n]->Clone();
  h_SMaxMjj[n]->Scale(sample[n].weight);

}


 TH2F *htmpMeta = new TH2F("htmpMeta","",40, 0, 8, 100, 0., 0.12);
 htmpMeta->SetXTitle("MaxEtajj");
 TH2F *htmpMjj = new TH2F("htmpMjj","",40, 0, 1000, 100, 0., 0.6); 
 htmpMjj->SetXTitle("MaxMjj");

 THStack *hs = new THStack("hs","Stacked 1D histograms");
 THStack *hs2 = new THStack("hs2","Stacked 1D histograms");
 TCanvas* c3 = new TCanvas("c3","",0,0,800,500);
 c3->Divide(2,1);
 c3->cd(1);
// htmpMeta->Draw();
 h_SMaxEta[0]->SetFillColor(kRed);
 h_SMaxEta[0]->SetFillStyle(3315);
 h_SMaxEta[0]->SetMarkerStyle(21);
 h_SMaxEta[2]->SetFillColor(kBlue);
 h_SMaxEta[2]->SetMarkerStyle(21); 
 hs->Add( h_SMaxEta[2]);
 hs->Add( h_SMaxEta[0]);
 hs->Draw();
 h_SMaxEta[1]->Draw("same");
 c3->cd(2);
 //htmpMjj->Draw();
 
 h_SMaxMjj[0]->SetFillColor(kRed);
 h_SMaxMjj[0]->SetFillStyle(3315);
 h_SMaxMjj[0]->SetMarkerStyle(21);
 h_SMaxMjj[2]->SetFillColor(kBlue);
 h_SMaxMjj[2]->SetMarkerStyle(21); 
 hs2->Add( h_SMaxMjj[2]);
 hs2->Add( h_SMaxMjj[0]);
 hs2->Draw();
 h_SMaxMjj[1]->Draw("same");
 
 c3->SaveAs("~/scratch0/c3.png");
/*

  TCanvas* c4 = new TCanvas("c4","",0,0,800,500);

 c4->Divide(2,1);
  c4->cd(1);
   htmpMeta->Draw();
  for(int n=0;n<nsam;n++)
  {

   if(n>runNMC) continue;
    h_MaxEta[n]->Scale(1/h_MaxEta[n]->GetEntries());
    h_MaxEta[n]->SetLineColor(n+4);
    h_MaxEta[n]->SetLineWidth(2);
    
    if(n==0){  h_MaxEta[n]->Draw("same");}
    else{ h_MaxEta[n]->Draw("same");}
   leg->AddEntry(h_MaxMjj[n],sample[n].tag , "pl"); 
  //cout<<"weight :"<< 1e9*weight[n]<<endl;
 // h_mcU2[n]->Write(); 
 // h_mcU2eta[n]->Write();

    

  }
  leg->Draw();
 
  c4->cd(2);
  htmpMjj->Draw();
   for(int n=0;n<nsam;n++)
   {
   if(n>runNMC) continue;
   h_MaxMjj[n]->Scale(1/h_MaxMjj[n]->GetEntries());
   h_MaxMjj[n]->SetLineColor(n+4);
   h_MaxMjj[n]->SetLineWidth(2); 
   if(n==0){  h_MaxMjj[n]->Draw("same");}
   else{ h_MaxMjj[n]->Draw("same");}
   }
  */



}
