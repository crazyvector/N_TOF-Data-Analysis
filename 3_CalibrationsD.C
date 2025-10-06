#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <set>
#include <algorithm>
#include <map>
#include <utility>

#include <TCanvas.h>
#include <TGraph.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TFile.h>
#include <TAxis.h>
#include <TFitResultPtr.h>

// Funcție ajutătoare pentru a salva pe o pagină nouă în PDF
void SaveCanvasToPDF(TCanvas *canvas, TString pdfFileName, bool &firstPage)
{
    if (firstPage)
    {
        canvas->SaveAs(pdfFileName + "(");
        firstPage = false;
    }
    else
    {
        canvas->SaveAs(pdfFileName);
    }
}

void fit(const char *fname, TString FitType, TString outputFolder)
{
    std::ifstream in(fname);
    if (!in.is_open())
    {
        std::cerr << "Error: Cannot open file " << fname << std::endl;
        return;
    }

    // Am creat o hartă pentru a stoca datele grupate pe Detector_ID.
    // Cheia este Detector_ID, iar valoarea este un vector de perechi (channel, energy).
    std::map<int, std::vector<std::pair<double, double>>> data_by_detector;
    std::string line;
    std::string source_name;
    int detector_id;
    double channel, energy, resolution;

    while (std::getline(in, line))
    {
        std::istringstream iss(line);
        // Datele din fișier sunt acum "Source DetectorID Channel Energy Resolution"
        if (iss >> source_name >> detector_id >> channel >> energy >> resolution)
        {
            data_by_detector[detector_id].push_back({channel, energy});
        }
    }
    in.close();

    // Verifică dacă am citit date din fișier.
    if (data_by_detector.empty())
    {
        std::cerr << "Error: No data found in file " << fname << std::endl;
        return;
    }

    TString filename = (FitType == "linear") ? "linear_fits" : "polynomial_fits";

    // Creez un singur canvas și un fișier PDF pentru toate rezultatele
    TCanvas *canvas = new TCanvas("canvas", filename, 800, 600);
    TString pdfFileName = outputFolder + "calibration_" + filename + ".pdf";
    bool firstPage = true;

    // Creez un singur fișier ROOT pentru toate rezultatele
    TFile *outputFile = new TFile(Form("%sEnergy_vs_Channels_All_Detectors_%s.root", outputFolder.Data(), FitType.Data()), "RECREATE");
    if (!outputFile || outputFile->IsZombie())
    {
        std::cerr << "Error: Cannot create output ROOT file." << std::endl;
        delete canvas;
        return;
    }
    outputFile->cd();

    // Iterează prin fiecare detector din hartă.
    for (auto const &[Detector_ID, data_points] : data_by_detector)
    {
        // Extrage datele pentru detectorul curent.
        int n = data_points.size();
        if (n < 2 && FitType == "linear")
        {
            std::cerr << "Error: Not enough points for fitting for Detector " << Detector_ID << ". Skipping.\n";
            continue;
        }
        if (n < 3 && FitType == "polynomial")
        {
            std::cerr << "Error: Not enough points for polynomial fit for Detector " << Detector_ID << ". Skipping.\n";
            continue;
        }

        std::vector<double> x, y;
        for (const auto &point : data_points)
        {
            x.push_back(point.first);
            y.push_back(point.second);
        }

        TString fitFormula;
        if (FitType == "linear")
            fitFormula = "[0] + [1]*x";
        else if (FitType == "polynomial")
            fitFormula = "[0] + [1]*x + [2]*x^2";

        // Grafic si functie
        TGraph *graph = new TGraph(n, &x[0], &y[0]);
        TF1 *fitFunc = new TF1("fitFunc", fitFormula, x[0], x[n - 1]);

        // Personalizare grafic
        graph->SetLineColor(kRed + 1);
        graph->SetLineWidth(3);
        graph->SetMarkerColor(kBlue + 2);
        graph->SetMarkerSize(1.5);
        graph->SetMarkerStyle(21);
        graph->SetTitle(Form("Energy vs Channels_Detector_%d", Detector_ID));
        graph->GetXaxis()->SetTitle("Channels");
        graph->GetYaxis()->SetTitle("Energy (keV)");

        // Fit si verificare
        TFitResultPtr fit = graph->Fit(fitFunc, "S");

        // Parametrii
        double intercept = fitFunc->GetParameter(0);                            // a
        double slope = fitFunc->GetParameter(1);                                // b
        double quad = (FitType == "polynomial") ? fitFunc->GetParameter(2) : 0; // c
        double chi2ndf = (fit->Ndf() > 0) ? fit->Chi2() / fit->Ndf() : -1;

        canvas->cd();
        graph->Draw("AP");
        fitFunc->Draw("SAME");

        TLegend *legend = new TLegend(0.15, 0.75, 0.55, 0.88);
        legend->AddEntry(graph, Form("Detector %d", Detector_ID), "p");
        legend->AddEntry((TObject *)0, Form("Slope = %.5f", slope), "");
        legend->AddEntry((TObject *)0, Form("Intercept = %.5f", intercept), "");
        legend->AddEntry((TObject *)0, Form("Quadratic = %.5f", quad), "");
        legend->AddEntry((TObject *)0, Form("#chi^{2}/NDF = %.4f", chi2ndf), "");
        legend->SetTextSize(0.03);
        legend->Draw();

        // Salvăm canvasul pe o pagină nouă în PDF
        SaveCanvasToPDF(canvas, pdfFileName, firstPage);

        // === Verificăm dacă există deja acest ID ===
        bool already_exists = false;
        std::ifstream in_check(Form("%scalibration_parameters.txt", outputFolder.Data()));
        if (in_check.is_open())
        {
            int existing_id;
            double s, i, d;
            std::string fit_type;
            while (in_check >> existing_id >> s >> i >> d >> fit_type)
            {
                if (existing_id == Detector_ID && fit_type == FitType)
                {
                    already_exists = true;
                    break;
                }
            }
            in_check.close();
        }

        // === Dacă nu există, scriem linia ===
        if (!already_exists)
        {
            std::ofstream out(Form("%scalibration_parameters.txt", outputFolder.Data()), std::ios::app);
            out << Detector_ID << " " << slope << " " << intercept << " " << quad << " " << FitType << std::endl;
            out.close();
        }
        else
        {
            std::cout << "⚠️ Detector " << Detector_ID << " already exists with fit type " << FitType << ". Skipping.\n";
        }

        // Salvăm graficul, funcția de fit și canvas-ul în fișierul ROOT principal
        graph->Write(Form("graph_Detector_%d", Detector_ID));
        fitFunc->Write(Form("fit_Detector_%d", Detector_ID));
        canvas->Write(Form("canvas_Detector_%d", Detector_ID));

        // Eliberarea memoriei pentru fiecare grafic și funcție
        delete graph;
        delete fitFunc;
    }

    // Închidem fișierul PDF
    canvas->Clear();
    canvas->SaveAs(pdfFileName + ")");
    delete canvas;

    // Închidem fișierul ROOT la finalul scriptului
    outputFile->Close();
    delete outputFile;
}
