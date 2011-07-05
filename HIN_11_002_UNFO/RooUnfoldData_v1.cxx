//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldData_v1.cxx,v 1.2 2011/04/21 12:44:08 yunju Exp $
//
// Description:
//      Simple example usage of the RooUnfold package using toy MC.
//
// Authors: Tim Adye <T.J.Adye@rl.ac.uk> and Fergus Wilson <fwilson@slac.stanford.edu>
//
//==============================================================================

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <iostream>
#include <stdio.h>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

//#include "RooUnfoldSvd.h"
//#include "RooUnfoldBinByBin.h"

#endif

//==============================================================================
// Unfolding HIN-11-002 Yun Ju Lu 2011 04 20
//==============================================================================



void RooUnfoldData(int setcbin=1,int method=1,std::string filename="HiFakeData.root",std::string hisname="measdata" )
{
#ifdef __CINT__
  gSystem->Load("../libRooUnfold.so");
#endif
gROOT->SetStyle("Plain");
gStyle->SetOptFit(0);
gStyle->SetOptDate(0);
gStyle->SetOptStat(0);
gStyle->SetPadRightMargin(0.05);
gStyle->SetPadTopMargin(0.08);
gStyle->SetStripDecimals(kTRUE);
gStyle->SetTickLength(0.03, "XYZ");
gStyle->SetPadTickX(1);
gStyle->SetPadTickY(1);
 

cout<<" ------------    Unfolding HIN-11-002 Yun Ju Lu 2011 04 20 ----------           "<<endl;
cout<<" ==============================================================================="<<endl;


// Read data file

   TH1D* hdata;
   TFile *f_data = TFile::Open(filename.data());
   hdata =(TH1D*)(f_data->Get(hisname.data()));




char chmet1[100]; 
if(method==1)
{

  sprintf(chmet1,"Bayes unfo");
} 
if(method==2)
{
   sprintf(chmet1,"Svd unfo ");
}
if(method==3)
{
 sprintf(chmet1,"BinByBin unfo");
}



printf(" Centrality bin : %d ; Method : %s \n",setcbin,chmet1);


  cout << "==================================== TRAIN ====================================" << endl;
    

     const int Nbins_P =10;
     double xbins[Nbins_P] = {5,10,15,20,25,30,40,50,80,140};
     //TCanvas *c11 = new TCanvas("c12","c12",600,500);
     //c11->cd(); 
     TH1D * tmp = (TH1D*)hdata->Clone();
     tmp->Reset();
     //tmp->Draw();
     //c11->Update();
     //cin.get();
     TH1D* hResTrue=(TH1D*)tmp->Clone() ;//Gen
     TH1D* hResMeas=(TH1D*)tmp->Clone() ; //Reco



    RooUnfoldResponse response (hResMeas,hResTrue);



    
     // ***************** histogram for testing sample *****************
    

     TH1D* hTrue= new TH1D ("true", "Test Truth",    Nbins_P-1, xbins); //Gen
     TH1D* hMeas= (TH1D*)tmp->Clone() ; //Reco

     //set centrality
     int Binmax;int Binmin;
     if(setcbin==1) {Binmax = 4 ;Binmin=0; }
     if(setcbin==2) {Binmax = 12 ;Binmin=4; }
     if(setcbin==3) {Binmax = 40 ;Binmin=12; }
    
    //******************* scale factor for each file ***************
    int scale15 = 15740 ;
    int scale30 = 1509;
    int scale50 = 240;
    int scale80 = 36; 

     //Pthat 15
     cout<<" Read Pt hat 15 - 30 "<<endl;

     
     TFile *ttfile = new TFile("HiRecoPhoton_pthat_15_30.root");
     TTree *ntuple =(TTree*)ttfile->Get("ntuple");
     
      double pt15respt[2];
     int pt15TcBin;

    TBranch *b_respt = ntuple->GetBranch("respt");
     b_respt->SetAddress(pt15respt);
     
     TBranch *b_TcBin = ntuple->GetBranch("TcBin");
     b_TcBin->SetAddress(&pt15TcBin); 
     



     Long64_t nentries = ntuple ->GetEntries();
     Long64_t nbytes = 0, nb = 0;



     for (Long64_t jentry=0; jentry<nentries;jentry++)
     {
        b_respt->GetEntry(jentry);
          b_TcBin->GetEntry(jentry);  
        if(pt15TcBin<Binmax&&pt15TcBin>=Binmin )
        {
           for(int nsca=0;nsca<scale15;nsca++)
              {  
                  response.Fill (pt15respt[1],pt15respt[0] );
              }
           /*
           if(jentry%2==1)
           {  
              for(int nsca=0;nsca<scale15;nsca++)
              {  
                  response.Fill (respt[1],respt[0] );
              }

           }
           else
           {
               for(int nsca=0;nsca<scale15;nsca++)
               {
               //hMeas->Fill(respt[1]);
               hTrue->Fill(respt[0]);
               }
           }
           */
       }

     }
     ttfile->Close(); 
    
     cout<<" Read Pt hat 30 - 50 "<<endl;

     //Pthat 30
     TFile *tfile2 = new TFile("HiRecoPhoton_pthat_30_50.root");
     TTree *ntuple2 =(TTree*)tfile2->Get("ntuple");
   
     double pt30respt[2];
     int pt30TcBin;
 
     TBranch *b30_respt = ntuple2->GetBranch("respt");
      b30_respt->SetAddress(pt30respt);

      TBranch *b30_TcBin = ntuple2->GetBranch("TcBin");
      b30_TcBin->SetAddress(&pt30TcBin);


     Long64_t nentries2 = ntuple2 ->GetEntries();
     Long64_t nbytes2 = 0, nb2 = 0;

  

     for (Long64_t jentry=0; jentry<nentries2;jentry++)
     {
      b30_respt->GetEntry(jentry);
        b30_TcBin->GetEntry(jentry); 

      
        if(pt30TcBin<Binmax&&pt30TcBin>=Binmin )
        {
           for(int nsca=0;nsca<scale30;nsca++)
              {  
                  response.Fill (pt30respt[1],pt30respt[0] );
              }
/* 
           if(jentry%2==1)
           {  
              for(int nsca=0;nsca<scale30;nsca++)
              {  
                  response.Fill (respt[1],respt[0] );
              }

           }
           else
           {
               for(int nsca=0;nsca<scale30;nsca++)
               {
              // hMeas->Fill(respt[1]);
               hTrue->Fill(respt[0]);
               }
           }
*/
       }

     }
   tfile2->Close(); 

     cout<<" Read Pt hat  50 - 80 "<<endl;

     //Pthat 50
     TFile *tfile3 = new TFile("HiRecoPhoton_pthat_50_80.root");
     TTree *ntuple3 =(TTree*)tfile3->Get("ntuple");
   
      double pt50respt[2];
      int pt50TcBin;

  //setting the branch address
  TBranch *b50_respt = ntuple3->GetBranch("respt");
  b50_respt->SetAddress(pt50respt);

  TBranch *b50_TcBin = ntuple3->GetBranch("TcBin");
  b50_TcBin->SetAddress(&pt50TcBin);
   

   


     Long64_t nentries3 = ntuple3 ->GetEntries();
     Long64_t nbytes3 = 0, nb3 = 0;

    
     for (Long64_t jentry=0; jentry<nentries3;jentry++)
     {
         b50_respt->GetEntry(jentry);
         b50_TcBin->GetEntry(jentry);
 

        if(pt50TcBin<Binmax&&pt50TcBin>=Binmin )
        {
           for(int nsca=0;nsca<scale50;nsca++)
           {  
                  response.Fill (pt50respt[1],pt50respt[0] );
           }
           /*
           if(jentry%2==1)
           {  
              for(int nsca=0;nsca<scale50;nsca++)
              {  
                  response.Fill (respt[1],respt[0] );
              }

           }
           else
           {
               for(int nsca=0;nsca<scale50;nsca++)
               {
              // hMeas->Fill(respt[1]);
               hTrue->Fill(respt[0]);
               }
           }
        */
       }

     }
   tfile3->Close(); 

     cout<<" Read Pt hat  80 - inf "<<endl;

     //Pthat 50
     TFile *tfile4 = new TFile("HiRecoPhoton_pthat_80_10000.root");
     TTree *ntuple4 =(TTree*)tfile4->Get("ntuple");
   
     double pt80respt[2];
      int pt80TcBin;

  //setting the branch address
    TBranch *b80_respt = ntuple4->GetBranch("respt");
    b80_respt->SetAddress(pt80respt);

    TBranch *b80_TcBin = ntuple4->GetBranch("TcBin");
    b80_TcBin->SetAddress(&pt80TcBin);
   

    
     Long64_t nentries4 = ntuple4 ->GetEntries();
     Long64_t nbytes4 = 0, nb4 = 0;

     for (Long64_t jentry=0; jentry<nentries4;jentry++)
     {
         b80_respt->GetEntry(jentry);
         b80_TcBin->GetEntry(jentry);
         
   

        if(pt80TcBin<Binmax&&pt80TcBin>=Binmin )
        {

           for(int nsca=0;nsca<scale80;nsca++)
           {  
             response.Fill (pt80respt[1],pt80respt[0] );
           }
           /*
           if(jentry%2==1)
           {  
              for(int nsca=0;nsca<scale80;nsca++)
              {  
                  response.Fill (respt[1],respt[0] );
              }

           }
           else
           {
               for(int nsca=0;nsca<scale80;nsca++)
               {
              // hMeas->Fill(respt[1]);
               hTrue->Fill(respt[0]);
               }
           }
           */

       }

     }
   tfile4->Close(); 





  cout << "==================================== TEST =====================================" << endl;


  cout << "==================================== UNFOLD ===================================" << endl;

char chmet[100]; 
if(method==1)
{
  RooUnfoldBayes    unfold (&response, hdata, 4);  
  sprintf(chmet,"Bayes unfo");
}  // OR
if(method==2)
{

RooUnfoldSvd      unfold (&response, hdata, 4);  
   sprintf(chmet,"Svd unfo ");
} // OR
if(method==3)
{

RooUnfoldBinByBin unfold (&response, hdata);
 sprintf(chmet,"BinByBin unfo");
}

   TH1D* hReco= (TH1D*) unfold.Hreco();
  // unfold.PrintTable (cout, hTrue);

   TCanvas *c10 = new TCanvas("c11","c11",600,500);
   TLegend *tleg = new TLegend(0.65, 0.6, 0.85, 0.85);
   tleg->SetBorderSize(0);
   TLatex *la1 = new TLatex(120, hReco->GetMaximum()*1.55, Form("Centrality %d",setcbin));
   tleg->AddEntry(hReco,chmet,"pl");
//   tleg->AddEntry(hTrue," gen ","pl");
   tleg->AddEntry(hdata," DATA ","pl");
   TH2F *htmp3 = new TH2F("htmp3","",150, 0., 150., 100, 0., hReco->GetMaximum()*1.5);
   htmp3->SetXTitle("E_{T} (GeV)");
   htmp3->Draw();
   hReco->SetMarkerColor(2);
   hReco->SetMarkerStyle(21);
   hReco->SetLineColor(2);
   hReco->Draw("SAME");
   hdata->Draw("histSAME");
  // hTrue->SetLineColor(8);
   //hTrue->Draw("histSAME");
   tleg->Draw("same");
   la1->Draw();
   c10->Print(Form("HiunfoPLots/centr_%d_Method_%d_v4.png",setcbin,method));
    printf("Result -----------  Centrality bin : %d ; Method : %s ------------ \n",setcbin,chmet1);
   for(int nbhis=0;nbhis<hdata->GetNbinsX();nbhis++)
   {
    
     printf("Bin center %5.2f  Input : %10.2f  ; Output %10.2f \n",hdata->GetBinCenter(nbhis+1),hdata->GetBinContent(nbhis+1),hReco->GetBinContent(nbhis+1)  );
   }
}


#ifndef __CINT__
int main () { RooUnfoldPhoJ(); return 0; }  // Main program when run stand-alone
#endif


