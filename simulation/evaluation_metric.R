evaluation_metric <-
  function(Net.true, Net.pred){
    
    K = dim(Net.true)[1]
    G = dim(Net.true)[2]
    
    
    TP = 0;
    FP = 0;
    TN = 0;
    FN = 0;
    
    for(k in 1:K){
      for(g in 1:G){
        Net.pred[[k,g]] = (abs(Net.pred[[k,g]]) > 1e-5)*upper.tri(Net.pred[[k,g]])
        Net.true[[k,g]] = Net.true[[k,g]]*upper.tri(Net.true[[k,g]])
        
        TP = TP + sum((Net.true[[k,g]]!=0)*(Net.pred[[k,g]]!=0))
        FP = FP + sum((Net.true[[k,g]]==0)*(Net.pred[[k,g]]!=0))
        TN = TN + sum((Net.true[[k,g]]==0)*(Net.pred[[k,g]]==0))
        FN = FN + sum((Net.true[[k,g]]!=0)*(Net.pred[[k,g]]==0))
        
      }
    }
    
    
    Pre = TP/(TP + FP+ 1e-6);
    Rec = TP / (TP + FN+1e-6);
    
    TPR = TP / (TP + FN+1e-6);
    FPR = FP / (FP + TN+1e-6);
    
    Per = c(Rec, Pre, FPR, TPR)
    
    result = list(Per=Per, Rec=Rec, Pre = Pre, FPR = FPR,  TPR= TPR)
    
  }