#include "TF1.h"

// histogram
int count(int Nev, int step, int SN, int ch, int Vin, int Vth)
{
   char* root = Form("../root/SN_%d/step_-%d/Vth_%d_Vin_%d_bd%d_ch0_Nevents_%d.root",SN, step, Vth, Vin, SN, Nev);
   //char* root = Form("step_-%d/Vth_%d_Vin_%d_bd%d_ch0_Nevents_%d.root",step, Vth, Vin, SN, Nev);

   TFile* f1 = new TFile(root);
   TTree* t = (TTree*)f1->Get("t");
   int nevents = t->GetEntries(Form("Sum$(tdc%d!=0)", ch));
   return nevents;
}

//fit function
//double fitf(int Vth, double mu, double sigma, int Nev)
//{
//  double Nhit;
//  Nhit = Nev*1.0 / (1 + exp((mu - Vth)/sigma));
//  return Nhit;
//}
double fitf(double* x, double* par)
{
   double mu = par[0];
   double sigma = par[1];
   double Nev = par[2];
   return  Nev / (1 + exp((mu - x[0])/sigma));
}

TCanvas* c1 = NULL;
void plot(int Nev, int step, int SN, int ch, int Vin, int Vth_min, int Vth_max)
{
   printf("Nev %d\n", Nev);
   printf("step %d\n", step);
   printf("SN %d\n", SN);
   printf("ch %d\n", ch);
   printf("Vin %d\n", Vin);
   printf("Vth_min %d\n", Vth_min);
   printf("Vth_max %d\n", Vth_max);

   gStyle->SetOptFit(1);

   int n=(Vth_max - Vth_min)/step + 1;
   double x[1000];
   double y[1000];
   double ex[1000];
   double ey[1000];
   double half_x;
   for (int i=0; i<n; i++) 
   {
      x[i] = Vth_min+i*step;
      y[i] = count(Nev, step, SN, ch, Vin, x[i]);
      ex[i] = 0.0;
      ey[i] = TMath::Sqrt(y[i]);
      if (y[i] > Nev/2*0.6 && y[i] < Nev/2*1.5) {
         half_x = x[i];
      }
   }
   printf("ch %d half_x %d\n", ch, half_x);
   TGraphErrors* gr = new TGraphErrors(n, x, y, ex, ey);
   gr->SetTitle(Form("Vth scan @ SN %d ch %d Vin %d; Vth; counts", SN, ch, Vin));
   gr->SetMarkerStyle(8);

   TF1*f = new TF1("f",fitf, Vth_min, Vth_max, 3);

   // center value
   f->SetParameter(0,half_x);
   f->SetParameter(1,1.5);
   f->FixParameter(2,Nev);

   if (c1==NULL) {
      c1 = new TCanvas("c1");
   }

   gr->Draw("ape");
   gr->Fit("f");
   c1->Modified();
   c1->Update();

   TPaveStats* stats = (TPaveStats*)c1->GetPrimitive("stats");
   stats->SetX1NDC(0.5);
   stats->SetX2NDC(0.86);
   stats->SetY1NDC(0.5);
   stats->SetY2NDC(0.7);
   stats->Draw();

   // save plot
   c1->Print(Form("../pdf/SN_%d/ch_%02d_Vin_%d.pdf", SN, ch, Vin));
   //c1->Print(Form("pdf/SN_%d/ch_%02d_Vin_%d.pdf", SN, ch, Vin));

   // write fit result
   FILE* fp = fopen(Form("../txt/SN_%d/ch_%02d_Vin_%d.txt", SN, ch, Vin), "w");
   //FILE* fp = fopen(Form("txt/SN_%d/ch_%02d_Vin_%d.txt", SN, ch, Vin), "w");
   double par[2];
   double err[2];
   for (int ipar=0; ipar<2; ipar++) {
      par[ipar] = f->GetParameter(ipar);
      err[ipar] = f->GetParError(ipar);
   }
   fprintf(fp, "SN %d ch %02d Vin %d Vth %f Vth_err %f sigma %f sigma_err %f\n", SN, ch, Vin, par[0], err[0], par[1], err[1]);
   fclose(fp);
}
void plotallch(int Nev, int step, int SN, int Vin, int Vth_min, int Vth_max)
{
   for (int ch=0; ch<48; ch++) {
      plot(Nev, step, SN, ch, Vin, Vth_min, Vth_max);
   }
}

void read_pars(int SN, int ch, double Vin, 
      double& Vth_read,
      double& Vth_err_read,
      double& sigma_read,
      double& sigma_err_read)
{
   int Vin_int = Vin;
   char line[1280];
   char* txt = Form("../txt/SN_%d/ch_%02d_Vin_%d.txt", SN, ch, Vin_int);
   //char* txt = Form("txt/SN_%d/ch_%02d_Vin_%d.txt", SN, ch, Vin_int);
   FILE* fp = fopen(txt, "r");
   if (fp==NULL) {
      fprintf(stderr,"cannot open %s\n", txt);
      return;
   }
   fgets(line, sizeof(line), fp);
   fclose(fp);
   printf("%s", line);

   int SN_read;
   int ch_read;
   int Vin_read;
   sscanf(line, "SN %d ch %02d Vin %d Vth %lf Vth_err %lf sigma %lf sigma_err %lf", 
         &SN_read, 
         &ch_read, 
         &Vin_read,
         &Vth_read,
         &Vth_err_read,
         &sigma_read,
         &sigma_err_read);
   printf("SN %d ch %02d Vin %d Vth %f Vth_err %f sigma %f sigma_err %f\n", 
         SN_read, 
         ch_read, 
         Vin_read, 
         Vth_read,
         Vth_err_read,
         sigma_read,
         sigma_err_read);
}
void fitch(int SN, int ch, int percent, double& offset, double& slope, double& offsetErr, double& slopeErr)
{
   gStyle->SetOptFit();

   double Vin[5] = {10, 20, 30, 40, 50};
   double VinErr[5] = {0};
   double parVth[5];
   double parSigma[5];
   double parVthErr[5];
   double parSigmaErr[5];

   for (int i=0; i<5; i++) {
      read_pars(SN, ch, Vin[i], parVth[i], parVthErr[i], parSigma[i], parSigmaErr[i]);
   }

   //for (int i=0; i<5; i++) {
   //   parVth[i] = 3800 - i*10;
   //   parVthErr[i] = 10;
   //}

   TGraphErrors* gr = new TGraphErrors(5, Vin, parVth, VinErr, parVthErr);
   gr->SetTitle(Form("Vth vs Vin @ SN %d ch %d; Vin (mV); Vth (mV)", SN, ch));
   gr->Draw("apl");

   gr->Fit("pol1");
   TF1*f1 = gr->GetFunction("pol1");
   offset = f1->GetParameter(0);
   slope  = f1->GetParameter(1);
   offsetErr = f1->GetParError(0);
   slopeErr  = f1->GetParError(1);
   printf("SN %d ch %02d offset %f slope %f offsetErr %f slopeErr %f\n", SN, ch, offset, slope, offsetErr, slopeErr);
}
void fitallch(int SN, int percent)
{
   double chs[48];
   double chsErr[48] = {0};
   double offset[48];
   double slope [48];
   double offsetErr[48];
   double slopeErr [48];
   char* hname1 = Form("SN%d_p%d_Offset", SN, percent);
   char* hname2 = Form("SN%d_p%d_Slope",  SN, percent);
   TH1F* h1 = new TH1F(hname1, hname1, 100, 3500, 4000);
   TH1F* h2 = new TH1F(hname2, hname2, 100, -2, 0);

   for (int ch=0; ch<48; ch++) {
      chs[ch] = ch;
      fitch(SN, ch, percent, offset[ch], slope[ch], offsetErr[ch], slopeErr[ch]);
      h1->Fill(offset[ch]);
      h2->Fill(slope[ch]);
   }
   TGraphErrors* gr1 = new TGraphErrors(48, chs, offset, chsErr, offsetErr);
   gr1->SetTitle(Form("Offset vs ch @ SN %d; ch; Offset (mV)", SN, ch));
   gr1->GetXaxis()->SetRangeUser(-1, 48);

   TGraphErrors* gr2 = new TGraphErrors(48, chs, slope, chsErr, slopeErr);
   gr2->SetTitle(Form("Slope vs ch @ SN %d; ch; Slope", SN, ch));
   gr2->GetXaxis()->SetRangeUser(-1, 48);

   TCanvas* c1 = new TCanvas("c1");
   c1->Divide(1,2);
   c1->cd(1); gr1->Draw("apl");
   c1->cd(2); gr2->Draw("apl");

   TCanvas* c2 = new TCanvas("c2");
   c2->Divide(1,2);
   c2->cd(1); h1->Draw();
   c2->cd(2); h2->Draw();

   c1->Print(Form("../fitallch/SN_%d/c1_%dpct.pdf",SN,percent));
   c2->Print(Form("../fitallch/SN_%d/c2_%dpct.pdf",SN,percent));
   //c1->Print(Form("fitallch/SN_%d/c1_%dpct.pdf",SN,percent));
   //c2->Print(Form("fitallch/SN_%d/c2_%dpct.pdf",SN,percent));
}
