## Module scores

# List the genes for the notochord module score
Noto_module1<- list(c("TBXT", "SHH", "NOG", "CHRD", "SOX9", "FOXA2"))

# Add the noto module scores to the NCC object

  NCC <- AddModuleScore(
    object = NCC,
    features = Noto_module1,
    ctrl = 5,
    name = 'Noto_module1')
  
# List the genes for the matrisome module score
  
  Noto_module4<- list(c("SPP1", "COL2A1", "HAPLN1", "HAPLN3", "SHH", "SEMA3C", "ANXA2", "CD109"))
  
  
  # add the matrisome module scores to the NCC object
  
  NCC <- AddModuleScore(
    object = NCC,
    features = Noto_module4,
    ctrl = 5,
    name = 'Noto_module4')
  
 
  # add the noto module score for integrated object
  
  CFS.EU.SCT.CC.clustered <- AddModuleScore(
    object = CFS.EU.SCT.CC.clustered,
    features = Noto_module1,
    ctrl = 5,
    name = 'Noto_module1')
  
  # add the matrisome module score for integrated object
  CFS.EU.SCT.CC.clustered <- AddModuleScore(
    object = CFS.EU.SCT.CC.clustered,
    features = Noto_module4,
    ctrl = 5,
    name = 'Noto_module4')
  

  #create dotplots for NCC and integrated object
  DotPlot(object = NCC, features = c("Noto_modeul11", "Noto_module41"))
  
  DotPlot(object = CFS.EU.SCT.CC.clustered, features = c("Noto_modeul11", "Noto_module41"))



  
  
  
  
  
  
