void checkResults_MuonGeneralAnalyzer(const char* fname="muonAnalysis.root"){
  printf("INFO: checking results in file \"%s\"\n", fname);
  TFile f(fname, "READ");
  if(! f.IsOpen())
    return;
  
  TH1F* counterOLD = (TH1F*) f.Get("muonAnalysisOLD/evt_muo_counter");
  TH1F* counterNEW = (TH1F*) f.Get("muonAnalysisNEW/evt_muo_counter");
  if(!counterOLD) { fprintf(stderr, "ERROR: muonAnalysisOLD/evt_muo_counter not found\n"); f.Close(); return ; }
  if(!counterNEW) { fprintf(stderr, "ERROR: muonAnalysisNEW/evt_muo_counter not found\n"); f.Close(); return ; }
  
  printf("%8.8s | %8.8s | %8.8s | %8.8s\n", "", "EVENTS", "MUONS", "RATIO");
  printf("%8.8s | %8.0f | %8.0f | %8.2f\n", "OLD",counterOLD->GetBinContent(1),counterOLD->GetBinContent(2),counterOLD->GetBinContent(2)/counterOLD->GetBinContent(1));
  printf("%8.8s | %8.0f | %8.0f | %8.2f\n", "NEW",counterNEW->GetBinContent(1),counterNEW->GetBinContent(2),counterNEW->GetBinContent(2)/counterNEW->GetBinContent(1));

  f.Close();
}
