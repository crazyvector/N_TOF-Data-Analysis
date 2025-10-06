#include <TROOT.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TAxis.h>
#include <TMath.h>
#include <TLine.h>
#include <TMarker.h>
#include <TText.h>
#include <TLegend.h>
#include <TFitResultPtr.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>

// --- Modelul cu 6 parametri (nemodificat)
double efficiencyFitFunc6(double *x, double *p) {
    double energy = x[0];
    if (energy <= 0) return 0;

    double lnE = TMath::Log(energy);
    double y = p[0] + p[1] * lnE + p[2] * TMath::Power(lnE, 2);
    double ln_epsilon = y * (2.0 / TMath::Pi()) *
                        TMath::ATan(TMath::Exp(p[3] + p[4] * lnE + p[5] * TMath::Power(lnE, 3))) - 25;

    return TMath::Exp(ln_epsilon);
}

// --- Modelul cu 4 parametri (nemodificat)
double efficiencyFitFunc4(double *x, double *p) {
    double energy = x[0];
    if (energy <= 0) return 0;

    double lnE = TMath::Log(energy);
    return TMath::Exp(p[0] + p[1] * lnE + p[2] * lnE * lnE + p[3] * lnE * lnE * lnE);
}

// !!! NOU: Funcția primește Detectorul și Sursa ca argumente
void FitEfficiency(const std::string &Detector, const std::string &Source, int modelChoice, double EnergyInterpolation) {

    std::string folder = "efficiency/";

    std::string pdfOut = folder + "EfficiencyFit_Det=" + Detector + "_Source=" + Source + ".png";
    std::string dataFile = folder + "PeakAreas_Detector=" + Detector + "_Source=" + Source + ".txt";

    gROOT->Reset();

    std::ifstream in(dataFile);
    if (!in.is_open()) {
        std::cerr << "❌ Nu pot deschide fișierul: " << dataFile << std::endl;
        return;
    }

    std::vector<double> energies, efficiencies, effErr;
    std::string line;
    bool firstLine = true;

    // Citirea datelor din fișierul generat de 6_PeakIntegrator.C
    while (std::getline(in, line)) {
        if (firstLine) { firstLine = false; continue; } // Sărim peste antet

        std::stringstream ss(line);
        std::string item;
        double energy = 0, eff = 0, effErr_ = 0;

        // Coloana 3: Energy, Coloana 26: Efficiency(%), Coloana 27: EfficiencyErr(%)
        for (int i = 0; i <= 28; ++i) {
            if (!std::getline(ss, item, ',')) break;
            if (i == 3)  energy = std::stod(item); // Folosim stod pentru conversie
            if (i == 26) eff = std::stod(item);
            if (i == 27) effErr_ = std::stod(item);
        }

        if (energy <= 0 || eff <= 0) continue;
        energies.push_back(energy);
        efficiencies.push_back(eff / 100.0); // Convertim % în fracție
        effErr.push_back(effErr_ / 100.0);
    }
    in.close();

    if (energies.empty()) {
        std::cerr << "❌ Nu s-au găsit puncte valide în fișierul " << dataFile << std::endl;
        return;
    }

    // --- Alegerea modelului ---
    TF1 *effFunc = nullptr;
    std::string funcFormula;

    if (modelChoice == 4) {
        effFunc = new TF1("effFunc", efficiencyFitFunc4, 50, 2000, 4);
        effFunc->SetParNames("p0", "p1", "p2", "p3");
        effFunc->SetParameters(-10.0, 2.0, -0.1, 0.01);
        funcFormula = "ε(E) = exp(p0 + p1*ln(E) + p2*ln²(E) + p3*ln³(E))";
    } else if (modelChoice == 6) {
        effFunc = new TF1("effFunc", efficiencyFitFunc6, 50, 2000, 6);
        effFunc->SetParNames("p0", "p1", "p2", "p3", "p4", "p5");
        effFunc->SetParameters(-10.0, 2.0, -0.1, -4.0, 0.2, -0.002);
        funcFormula = "ln(ε) = y * (2/π) * atan(exp(...)) - 25";
    } else {
         std::cerr << "❌ Model de fit invalid. Se folosește modelul implicit (4 parametri)." << std::endl;
         effFunc = new TF1("effFunc", efficiencyFitFunc4, 50, 2000, 4);
         modelChoice = 4;
    }

    effFunc->SetLineColor(kRed);
    effFunc->SetLineWidth(2);

    // --- Graficul ---
    TGraphErrors *g = new TGraphErrors(energies.size(), &energies[0], &efficiencies[0], 0, &effErr[0]);
    g->SetTitle(Form("Efficiency Curve for %s (%s);Energy [keV];Efficiency", Source.c_str(), Detector.c_str()));
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kBlue);

    TCanvas *c1 = new TCanvas("c1", "Efficiency Fit", 800, 600);

    g->Draw("AP");

    // --- Fit-ul ---
    TFitResultPtr r = g->Fit(effFunc, "RS");

    // --- Legendă ---
    TLegend *legend = new TLegend(0.45, 0.75, 0.88, 0.88);
    legend->SetTextSize(0.025);
    legend->SetBorderSize(1);
    legend->SetFillColor(0);
    legend->AddEntry(g, "Data points", "p");

    if (modelChoice == 4) {
        legend->AddEntry(effFunc, "Eff(E) = exp(p0 + p1*ln(E)", "");
        legend->AddEntry((TObject*)0, " + p2*ln(E)^{2} + p3*ln(E)^{3})", "");}
    else {
        legend->AddEntry(effFunc, "ln(eps) = (p0 + p1*ln(E) + p2*ln(E)^{2}) * (2/#pi)", "");
        legend->AddEntry((TObject*)0, " * atan(exp(p3 + p4*ln(E) + p5*ln(E)^{3})) - 25", "");
    }

    legend->Draw();

    if (r->IsValid()) {
        std::cout << "✅ Fit reușit!\n📈 Parametrii obținuți:\n";
        for (int i = 0; i < effFunc->GetNpar(); ++i) {
            std::cout << Form("  p%d = %.5f +/- %.5f", i, effFunc->GetParameter(i), effFunc->GetParError(i)) << std::endl;
        }

        // --- Interpolare la 1368 keV ---
        double xVal = EnergyInterpolation;
        double yVal = effFunc->Eval(xVal);

        TLine *line = new TLine(xVal, g->GetYaxis()->GetXmin(), xVal, yVal);
        line->SetLineColor(kGreen + 2);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw();

        TMarker *marker = new TMarker(xVal, yVal, 21);
        marker->SetMarkerColor(kGreen + 2);
        marker->SetMarkerSize(1.2);
        marker->Draw();

        TPaveText *pt = new TPaveText(0.15, 0.75, 0.40, 0.88, "NDC");
        pt->SetFillStyle(0);
        pt->SetBorderSize(0);
        pt->AddText(Form("Interpolated at %.0f keV:", xVal));
        pt->AddText(Form("Eff = %.4f (%.4f %%)", yVal, yVal * 100));
        pt->Draw();

    } else {
        std::cerr << "❌ Fit eșuat." << std::endl;
    }

    c1->SaveAs(pdfOut.c_str());

    // Se salvează parametrii de fit într-un fișier, dacă este necesar (adăugare opțională)
    std::string name = folder + "efficiency_fit_params.txt";
    std::ofstream fitParamsFile(name);

    delete c1;
    delete g;
    delete effFunc;
}