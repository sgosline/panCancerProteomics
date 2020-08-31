##train and test model
require(e1071)
library(readxl)
require(dplyr)
train<-function(){
    all.dat<-readxl::read_xlsx('inst/41590_2017_BFni3693_MOESM10_ESM.xlsx')
    activated<-all.dat%>%
      dplyr::select("Majority protein IDs",'Gene names',starts_with('Intensity'))%>%
      dplyr::select("Gene names","Majority protein IDs",ends_with('activated'))%>%
      tidyr::pivot_longer(cols=ends_with('activated'),values_to='Activated',names_to='Sample')%>%
      rowwise()%>%mutate(Sample=stringr::str_replace(Sample,'_activated',''))%>%
      #mutate(geneSamp=paste(`Gene names`,Sample))%>%
     ungroup()%>%distinct()
    #  select(Activated,geneSamp)

    steady<-all.dat%>%
      dplyr::select("Majority protein IDs","Gene names",starts_with('Intensity'))%>%
      dplyr::select("Majority protein IDs","Gene names",ends_with('steady-state'))%>%
      tidyr::pivot_longer(cols=ends_with('steady-state'),values_to='Steady',names_to='Sample')%>%
      rowwise()%>%mutate(Sample=stringr::str_replace(Sample,'_steady-state',''))%>%
      #mutate(geneSamp=paste(`Gene names`,Sample))%>%
      ungroup()%>%
      distinct()
      #select(Steady,geneSamp)s

    combined<-activated%>%inner_join(steady)%>%
      mutate(logRatio=log10(0.01+Activated)-log10(0.01+Steady))%>%
      #rowwise()%>%mutate(Sample=stringr::str_replace(Sample,'Intensity_',''))%>%
      separate(Sample,into=c("data","cellType","replicate"),sep='_',remove=FALSE)

    dat<-combined%>%select(`Gene names`,Sample,logRatio,cellType)%>%
      pivot_wider(names_from='Gene names',values_from=logRatio,values_fn=list(logRatio=mean),
                  values_fill=list(logRatio=0.0))%>%
      tibble::column_to_rownames('Sample')
    # yvar<-combined%>%select(Sample,cellType)%>%distinct()
    zvar<-which(apply(dat,2,var)==0)
    if(length(zvar)>0)
      dat<-dat[,-zvar]
    mod=svm(cellType~.,data=dat)
  return(mod)
}
