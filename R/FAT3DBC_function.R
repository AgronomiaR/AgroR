#' Analysis: DBC experiments in triple factorial
#'
#' @description Analysis of an experiment conducted in a randomized block design in a triple factorial scheme using analysis of variance of fixed effects.
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @author Leandro Simoes Azeredo Goncalves
#' @author Rodrigo Yudi Palhaci Marubayashi
#' @param f1 Numeric or complex vector with factor 1 levels
#' @param f2 Numeric or complex vector with factor 2 levels
#' @param f3 Numeric or complex vector with factor 3 levels
#' @param block Numerical or complex vector with blocks
#' @param response Numerical vector containing the response of the experiment.
#' @param mcomp Multiple comparison test (Tukey (\emph{default}), LSD, Scott-Knott and Duncan)
#' @param quali Defines whether the factor is quantitative or qualitative (\emph{qualitative})
#' @param names.fat Allows labeling the factors 1, 2 and 3.
#' @param grau Degree of polynomial in case of quantitative factor (\emph{default} is 1)
#' @param xlab Treatments name (Accepts the \emph{expression}() function)
#' @param ylab Variable response name (Accepts the \emph{expression}() function)
#' @param alpha.t Significance level of the multiple comparison test (\emph{default} is 0.05)
#' @param alpha.f Level of significance of the F test (\emph{default} is 0.05)
#' @param norm Error normality test (\emph{default} is Shapiro-Wilk)
#' @param homog Homogeneity test of variances (\emph{default} is Bartlett)
#' @param transf Applies data transformation (\emph{default} is 1; for log consider 0)
#' @param sup Number of units above the standard deviation or average bar on the graph
#' @param geom Graph type (columns or segments)
#' @param fill Defines chart color (to generate different colors for different treatments, define fill = "trat")
#' @param angulo x-axis scale text rotation
#' @param textsize Font size
#' @param dec Number of cells
#' @param family Font family
#' @param theme ggplot2 theme (\emph{default} is theme_classic())
#' @param addmean Plot the average value on the graph (\emph{default} is TRUE)
#' @param errorbar Plot the standard deviation bar on the graph (In the case of a segment and column graph) - \emph{default} is TRUE
#' @param angle.label label angle
#' @return The analysis of variance table, the Shapiro-Wilk error normality test, the Bartlett homogeneity test of variances, the Durbin-Watson error independence test, multiple comparison test (Tukey, LSD, Scott-Knott or Duncan) or adjustment of regression models up to grade 3 polynomial, in the case of quantitative treatments. The column chart for qualitative treatments is also returned.For significant triple interaction only, no graph is returned.
#' @note The function does not perform multiple regression in the case of two or more quantitative factors. The bars of the column and segment graphs are standard deviation.
#' @note In the final output when transformation (transf argument) is different from 1, the columns resp and respo in the mean test are returned, indicating transformed and non-transformed mean, respectively.
#' @references
#'
#' Principles and procedures of statistics a biometrical approach Steel, Torry and Dickey. Third Edition 1997
#'
#' Multiple comparisons theory and methods. Departament of statistics the Ohio State University. USA, 1996. Jason C. Hsu. Chapman Hall/CRC.
#'
#' Practical Nonparametrics Statistics. W.J. Conover, 1999
#'
#' Ramalho M.A.P., Ferreira D.F., Oliveira A.C. 2000. Experimentacao em Genetica e Melhoramento de Plantas. Editora UFLA.
#'
#' Scott R.J., Knott M. 1974. A cluster analysis method for grouping mans in the analysis of variance. Biometrics, 30, 507-512.
#'
#' Ferreira, E. B., Cavalcanti, P. P., and Nogueira, D. A. (2014). ExpDes: an R package for ANOVA and experimental designs. Applied Mathematics, 5(19), 2952.
#'
#' Mendiburu, F., and de Mendiburu, M. F. (2019). Package ‘agricolae’. R Package, Version, 1-2.
#'
#' @keywords DIC
#' @keywords Factorial
#' @export
#' @examples
#' library(AgroR)
#' data(enxofre)
#' with(enxofre, FAT3DBC(f1, f2, f3, bloco, resp))

FAT3DBC=function(f1,
                 f2,
                 f3,
                 block,
                 response,
                 quali=c(TRUE,TRUE,TRUE),
                 names.fat=c("F1","F2","F3"),
                 mcomp='tukey',
                 transf=1,
                 norm="sw",
                 homog="bt",
                 alpha.f=0.05,
                 alpha.t=0.05,
                 ylab="Response",
                 xlab="",
                 sup=NA,
                 grau=NA,
                 fill="lightblue",
                 theme=theme_classic(),
                 angulo=0,
                 errorbar=TRUE,
                 addmean=TRUE,
                 family="sans",
                 dec=3,
                 geom="bar",
                 textsize=12,
                 angle.label=0) {
    if(is.na(sup==TRUE)){sup=0.2*mean(response)}
    if(angle.label==0){hjust=0.5}else{hjust=0}
    fator1=f1
    fator2=f2
    fator3=f3
    fator1a=fator1
    fator2a=fator2
    fator3a=fator3


    bloco=block
    fac.names=names.fat
    requireNamespace("ScottKnott")
    requireNamespace("crayon")
    requireNamespace("ggplot2")
    requireNamespace("nortest")
    # ================================
    # Transformacao de dados
    # ================================
    if(transf==1){resp=response}else{resp=(response^transf-1)/transf}
    if(transf==0){resp=log(response)}
    if(transf==0.5){resp=sqrt(response)}
    if(transf==-0.5){resp=1/sqrt(response)}
    if(transf==-1){resp=1/response}

    # =================================
    ## Reconstruindo vetores
    # =================================
    fatores<-data.frame(fator1,fator2,fator3)
    Fator1<-factor(fator1,levels=unique(fator1));
    Fator2<-factor(fator2,levels=unique(fator2));
    Fator3<-factor(fator3,levels=unique(fator3))
    nv1<-length(summary(Fator1)); nv2<-length(summary(Fator2)); nv3<-length(summary(Fator3))
    J<-(length(resp))/(nv1*nv2*nv3)
    lf1<-levels(Fator1); lf2<-levels(Fator2); lf3<-levels(Fator3)
    bloco=as.factor(bloco)

    # =================================
    ## Anova
    # =================================
    anava<-aov(resp~Fator1*Fator2*Fator3+bloco)
    anavaF3<-anova(anava)
    anovaF3=anavaF3
    colnames(anovaF3)=c("GL","SQ","QM","Fcal","p-value")
    respad=anava$residuals/sqrt(anavaF3$`Mean Sq`[9])
    out=respad[respad>3 | respad<(-3)]
    out=names(out)
    out=if(length(out)==0)("No discrepant point")else{out}
    resids=anava$residuals/sqrt(anavaF3$`Mean Sq`[9])
    Ids=ifelse(resids>3 | resids<(-3), "darkblue","black")
    residplot=ggplot(data=data.frame(resids,Ids),aes(y=resids,x=1:length(resids)))+
        geom_point(shape=21,color="gray",fill="gray",size=3)+
        labs(x="",y="Standardized residuals")+
        geom_text(x=1:length(resids),label=1:length(resids),color=Ids,size=4)+
        scale_x_continuous(breaks=1:length(resids))+
        theme_classic()+theme(axis.text.y = element_text(size=12),
                              axis.text.x = element_blank())+
        geom_hline(yintercept = c(0,-3,3),lty=c(1,2,2),color="red",size=1)

    # =================================
    ## Saída inicial
    # =================================

    #Teste de normalidade
    norm1<-shapiro.test(anava$residuals)
    cat(green(bold("\n------------------------------------------")))
    cat(green(bold("\nNormality of errors")))
    cat(green(bold("\n------------------------------------------")))
    print(norm1)
    message(if(norm1$p.value>0.05){
        black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, errors can be considered normal")}
        else {"As the calculated p-value is less than the 5% significance level, H0 is rejected. Therefore, errors do not follow a normal distribution"})

    # Teste de homogeneidade das variâncias
    homog1=bartlett.test(anava$residuals~paste(Fator1,Fator2,Fator3))
    cat(green(bold("\n------------------------------------------")))
    cat(green(bold("\nHomogeneity of Variances")))
    cat(green(bold("\n------------------------------------------")))
    print(homog1)
    message(if(homog1$p.value[1]>0.05){
        black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, the variances can be considered homogeneous")}
        else {"As the calculated p-value is less than the 5% significance level, H0 is rejected. Therefore, the variances are not homogeneous"})

    # Independencia dos erros
    indep=dwtest(anava)
    cat(green(bold("\n------------------------------------------")))
    cat(green(bold("\nIndependence from errors")))
    cat(green(bold("\n------------------------------------------")))
    print(indep)
    message(if(indep$p.value>0.05){
        black("As the calculated p-value is greater than the 5% significance level, hypothesis H0 is not rejected. Therefore, errors can be considered independent")}
        else {"As the calculated p-value is less than the 5% significance level, H0 is rejected. Therefore, errors are not independent"})
    cat("\n")
    cat(green(bold("\n-----------------------------------------------------------------\n")))
    cat(green(bold("Additional Information")))
    cat(green(bold("\n-----------------------------------------------------------------\n")))
    cat(paste("\nCV (%) = ",round(sqrt(anavaF3$`Mean Sq`[9])/mean(resp,na.rm=TRUE)*100,2)))
    cat(paste("\nMean = ",round(mean(response,na.rm=TRUE),4)))
    cat(paste("\nMedian = ",round(median(response,na.rm=TRUE),4)))
    cat("\nPossible outliers = ", out)
    cat("\n")
    cat(green(bold("\n------------------------------------------\n")))
    cat(green(bold("Analysis of Variance")))
    cat(green(bold("\n------------------------------------------\n")))
    anava1=as.matrix(data.frame(anovaF3))
    colnames(anava1)=c("Df","Sum Sq","Mean.Sq","F value","Pr(F)" )
    print(anava1,na.print = "")
    cat("\n")

    if(transf==1 && norm1$p.value<0.05 | transf==1 && indep$p.value<0.05 | transf==1 &&homog1$p.value<0.05){
        message("\n Your analysis is not valid, suggests using a non-parametric test and try to transform the data\n")}else{}
    if(transf != 1 && norm1$p.value<0.05 | transf!=1 && indep$p.value<0.05 | transf!=1 && homog1$p.value<0.05){
        message("\n Your analysis is not valid, suggests using the function FATDIC.art\n")}else{}
    message(if(transf !=1){blue("\nNOTE: resp = transformed means; respO = averages without transforming\n")})

    ################################################################################################
    # Efeitos simples
    ################################################################################################
    if(anavaF3[5,5]>alpha.f && anavaF3[6,5]>alpha.f && anavaF3[7,5]>alpha.f && anavaF3[8,5]>alpha.f) {
        graficos=list(1,2,3)
        cat(green(bold("\n------------------------------------------\n")))
        cat(green(bold('Non-significant interaction: analyzing the simple effects')))
        cat(green(bold("\n------------------------------------------\n")))
        fatores<-data.frame('fator 1'=fator1,'fator 2' = fator2,'fator 3' = fator3)

        for(i in 1:3){
            # Comparação múltipla
            if(quali[i]==TRUE && anavaF3[i,5]<=alpha.f) {
                cat(green(bold("\n------------------------------------------\n")))
                cat(fac.names[i])
                cat(green(bold("\n------------------------------------------\n")))
                if(mcomp=='tukey'){letra=TUKEY(resp,fatores[,i],anavaF3[9,1],anavaF3[9,3],alpha.t)
                letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
                if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
                if(mcomp=="sk"){
                    ad=data.frame(Fator1,Fator2,Fator3)
                    letra=SK(anava,colnames(ad[i]))
                    letra1=data.frame(resp=letra$m.inf[,1],groups=letters[letra$groups])
                    #letra1$resp=as.numeric(as.character(letra1$resp))
                    if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
                if(mcomp=="duncan"){
                    ad=data.frame(Fator1,Fator2,Fator3)
                    letra <- duncan(anava, colnames(ad[i]), alpha=alpha.t)
                    letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
                    if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
                if(mcomp=="lsd"){
                    ad=data.frame(Fator1,Fator2,Fator3)
                    letra <- LSD(anava, colnames(ad[i]), alpha=alpha.t)
                    letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
                    if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
                print(letra1)
                cat(green(bold("\n------------------------------------------\n")))
                dadosm=data.frame(letra1,
                                  media=tapply(response, c(fatores[i]), mean, na.rm=TRUE)[rownames(letra1)],
                                  desvio=tapply(response, c(fatores[i]), sd, na.rm=TRUE)[rownames(letra1)])
                dadosm$Tratamentos=factor(rownames(dadosm),levels = unique(unlist(fatores[i])))
                dadosm$limite=dadosm$media+dadosm$desvio
                dadosm=dadosm[as.character(unique(unlist(fatores[i]))),]
                if(addmean==TRUE){dadosm$letra=paste(format(dadosm$media,digits = dec),dadosm$groups)}
                if(addmean==FALSE){dadosm$letra=dadosm$groups}
                media=dadosm$media
                desvio=dadosm$desvio
                Tratamentos=dadosm$Tratamentos
                letra=dadosm$letra
                if(geom=="bar"){grafico=ggplot(dadosm,
                                               aes(x=Tratamentos,
                                                   y=media))
                if(fill=="trat"){grafico=grafico+
                    geom_col(aes(fill=Tratamentos),
                             color=1)}else{grafico=grafico+
                                 geom_col(aes(fill=Tratamentos),
                                          fill=fill,color=1)}
                if(errorbar==TRUE){grafico=grafico+
                    geom_text(aes(y=media+sup+
                                      if(sup<0){-desvio}else{desvio},
                                  label=letra),family=family,angle=angle.label, hjust=hjust)}
                if(errorbar==FALSE){grafico=grafico+
                    geom_text(aes(y=media+sup,label=letra),family=family,angle=angle.label, hjust=hjust)}
                if(errorbar==TRUE){grafico=grafico+
                    geom_errorbar(data=dadosm,
                                  aes(ymin=media-desvio,
                                      ymax=media+desvio,color=1),
                                  color="black",width=0.3)}
                grafico=grafico+theme+
                    ylab(ylab)+
                    xlab(xlab)+
                    theme(text = element_text(size=textsize,color="black", family = family),
                          axis.text = element_text(size=textsize,color="black", family = family),
                          axis.title = element_text(size=textsize,color="black", family = family),
                          legend.position = "none")
                print(grafico)}

                # ================================
                # grafico de segmentos
                # ================================
                if(geom=="point"){grafico=ggplot(dadosm,
                                                 aes(x=Tratamentos,
                                                     y=media))
                if(fill=="trat"){grafico=grafico+
                    geom_point(aes(color=Tratamentos))}
                else{grafico=grafico+
                    geom_point(aes(color=Tratamentos),color=fill,size=4)}
                if(errorbar==TRUE){grafico=grafico+
                    geom_text(aes(y=media+sup+if(sup<0){-desvio}else{desvio},
                                  label=letra),family=family,angle=angle.label, hjust=hjust)}
                if(errorbar==FALSE){grafico=grafico+
                    geom_text(aes(y=media+sup,
                                  label=letra),family=family,angle=angle.label, hjust=hjust)}
                if(errorbar==TRUE){grafico=grafico+
                    geom_errorbar(data=dadosm,
                                  aes(ymin=media-desvio,
                                      ymax=media+desvio,color=1),
                                  color="black",width=0.3)}
                grafico=grafico+theme+
                    ylab(ylab)+
                    xlab(xlab)+
                    theme(text = element_text(size=textsize,color="black", family = family),
                          axis.text = element_text(size=textsize,color="black", family = family),
                          axis.title = element_text(size=textsize,color="black", family = family),
                          legend.position = "none")
                print(grafico)}
            }

            if(anavaF3[i,5]>alpha.f) {
                cat(green(bold("\n------------------------------------------\n")))
                cat(fac.names[i])
                cat(green(bold("\n------------------------------------------\n")))
                mean.table<-mean.stat(response,fatores[,i],mean)
                colnames(mean.table)<-c('Niveis','Medias')
                print(mean.table)
                grafico=NA}

            # Regressão
            if(quali[i]==FALSE && anavaF3[i,5]<=alpha.f){
                cat(fac.names[i])
                dose=as.numeric(as.character(as.vector(unlist(fatores[,i]))))
                grafico=polynomial(dose,resp,grau = grau)
                cat(green("To edit graphical parameters, I suggest analyzing using the \"polynomial\" command"))
                cat(green(bold("\n------------------------------------------")))}
            graficos[[1]]=residplot
            graficos[[i+1]]=grafico
        }
    }

    #####################################################################
    #Interacao Fator1*Fator2      +     Fator3
    #####################################################################
    if(anavaF3[8,5]>alpha.f && anavaF3[5,5]<=alpha.f){
    cat(green(bold("\n------------------------------------------\n")))
        cat(green(bold("Interaction",paste(fac.names[1],'*',fac.names[2],sep='')," significant: unfolding the interaction")))
    cat(green(bold("\n------------------------------------------\n")))

        #Desdobramento de FATOR 1 dentro do niveis de FATOR 2
        cat(green(bold("\n------------------------------------------\n")))
        cat("Analyzing ", fac.names[1], ' within the combination of levels ', fac.names[2])
        cat(green(bold("\n------------------------------------------\n")))
        des<-aov(resp~Fator2/Fator1+Fator3+Fator2+Fator2:Fator3+Fator1:Fator2:Fator3+bloco)
        l<-vector('list',nv2)
        names(l)<-names(summary(Fator2))
        v<-numeric(0)
        for(j in 1:nv2) {
            for(i in 0:(nv1-2)) v<-cbind(v,i*nv2+j)
            l[[j]]<-v
            v<-numeric(0)}
        des1<-summary(des,split=list('Fator2:Fator1'=l))[[1]]
        des1a=des1[-c(1,2,3,length(des1[,1]),length(des1[,1])-1,length(des1[,1])-2),]
        print(des1a)

        # Teste de Tukey
        if(quali[1]==TRUE & quali[2]==TRUE){
            if (mcomp == "tukey"){
                tukeygrafico=c()
                ordem=c()
                for (i in 1:nv2) {
                    trati=fatores[, 1][Fator2 == lf2[i]]
                    trati=factor(trati,levels = unique(trati))
                    respi=resp[Fator2 == lf2[i]]
                    mod=aov(respi~trati)
                    tukey=TUKEY(mod,"trati",anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                    if(transf !="1"){tukey$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(tukey$groups)]}
                    tukeygrafico[[i]]=tukey$groups[levels(trati),2]
                    ordem[[i]]=rownames(tukey$groups[levels(trati),])
                    }
                letra=unlist(tukeygrafico)
                datag=data.frame(letra,ordem=unlist(ordem))
                datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
                datag=datag[order(datag$ordem),]
                letra=datag$letra}
            if (mcomp == "duncan"){
                duncangrafico=c()
                ordem=c()
                for (i in 1:nv2) {
                    trati=fatores[, 1][Fator2 == lf2[i]]
                    trati=factor(trati,levels = unique(trati))
                    respi=resp[Fator2 == lf2[i]]
                    duncan=duncan(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                    if(transf !="1"){duncan$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(duncan$groups)]}
                    duncangrafico[[i]]=duncan$groups[levels(trati),2]
                    ordem[[i]]=rownames(duncan$groups[levels(trati),])
                    }
                letra=unlist(duncangrafico)
                datag=data.frame(letra,ordem=unlist(ordem))
                datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
                datag=datag[order(datag$ordem),]
                letra=datag$letra}
            if (mcomp == "lsd"){
                lsdgrafico=c()
                ordem=c()
                for (i in 1:nv2) {
                    trati=fatores[, 1][Fator2 == lf2[i]]
                    trati=factor(trati,levels = unique(trati))
                    respi=resp[Fator2 == lf2[i]]
                    lsd=LSD(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                    if(transf !="1"){lsd$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(lsd$groups)]}
                    lsdgrafico[[i]]=lsd$groups[levels(trati),2]
                    ordem[[i]]=rownames(lsd$groups[levels(trati),])
                    }
                letra=unlist(lsdgrafico)
                datag=data.frame(letra,ordem=unlist(ordem))
                datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
                datag=datag[order(datag$ordem),]
                letra=datag$letra}
            if (mcomp == "sk"){
                skgrafico=c()
                ordem=c()
                for (i in 1:nv2) {
                    trati=fatores[, 1][Fator2 == lf2[i]]
                    trati=factor(trati,levels = unique(trati))
                    respi=resp[Fator2 == lf2[i]]
                    sk=sk(respi,trati,anavaF3$Df[9],anavaF3$`Sum Sq`[9],alpha.t)
                    if(transf !="1"){sk$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(sk$groups)]}
                    skgrafico[[i]]=sk[levels(trati),2]
                    ordem[[i]]=rownames(sk[levels(trati),])
                    }
                letra=unlist(skgrafico)
                datag=data.frame(letra,ordem=unlist(ordem))
                datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
                datag=datag[order(datag$ordem),]
                letra=datag$letra}}

        # Desdobramento de F2 dentro de F1

        cat(green(bold("\n------------------------------------------\n")))
        cat("Analyzing ", fac.names[2], " inside of the level of ",fac.names[1])
        cat(green(bold("\n------------------------------------------\n")))

        # Desdobramento de F1 dentro de F2
        des<-aov(resp~Fator1/Fator2+Fator3+Fator1+Fator1:Fator3+Fator1:Fator2:Fator3)
        l<-vector('list',nv1)
        names(l)<-names(summary(Fator1))
        v<-numeric(0)
        for(j in 1:nv1) {
                for(i in 0:(nv2-2)) v<-cbind(v,i*nv1+j)
                l[[j]]<-v
                v<-numeric(0)}
        des1<-summary(des,split=list('Fator1:Fator2'=l))[[1]]
        des1a=des1[-c(1,2,3,length(des1[,1]),length(des1[,1])-1,length(des1[,1])-2),]
        print(des1a)

        #-------------------------------------
        # Teste de Tukey
        #-------------------------------------
        if(quali[1]==TRUE & quali[2]==TRUE){
                if (mcomp == "tukey"){
                    tukeygrafico1=c()
                    for (i in 1:nv1) {
                        trati=as.factor(fatores[, 2][Fator1 == lf1[i]])
                        respi=resp[Fator1 == lf1[i]]
                        mod=aov(respi~trati)
                        tukey=TUKEY(mod,"trati",anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                        if(transf !="1"){tukey$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(tukey$groups)]}
                        tukeygrafico1[[i]]=tukey$groups[levels(trati),2]
                        }
                    letra1=unlist(tukeygrafico1)
                    letra1=toupper(letra1)}
                if (mcomp == "duncan"){
                    duncangrafico1=c()
                    for (i in 1:nv1) {
                        trati=as.factor(fatores[, 2][Fator1 == lf1[i]])
                        respi=resp[Fator1 == lf1[i]]
                        duncan=duncan(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                        if(transf !="1"){duncan$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(duncan$groups)]}
                        duncangrafico1[[i]]=duncan$groups[levels(trati),2]
                        }
                    letra1=unlist(duncangrafico1)
                    letra1=toupper(letra1)}
                if (mcomp == "lsd"){
                    lsdgrafico1=c()
                    for (i in 1:nv1) {
                        trati=as.factor(fatores[, 2][Fator1 == lf1[i]])
                        respi=resp[Fator1 == lf1[i]]
                        lsd=LSD(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                        if(transf !="1"){lsd$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(lsd$groups)]}
                        lsdgrafico1[[i]]=lsd$groups[levels(trati),2]
                        }
                    letra1=unlist(lsdgrafico1)
                    letra1=toupper(letra1)}
                if (mcomp == "sk"){
                    skgrafico1=c()
                    for (i in 1:nv1) {
                        trati=fatores[, 2][Fator1 == lf1[i]]
                        trati=factor(trati,levels = unique(trati))
                        respi=resp[Fator1 == lf1[i]]
                        sk=sk(respi,trati,anavaF3$Df[9],anavaF3$`Sum Sq`[9],alpha.t)
                        if(transf !="1"){sk$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(sk)]}
                        skgrafico1[[i]]=sk[levels(trati),2]
                        }
                    letra1=unlist(skgrafico1)
                    letra1=toupper(letra1)}}

        # -----------------------------
        # Gráfico de colunas
        #------------------------------
        if(quali[1] & quali[2]==TRUE){
                    f1=rep(levels(Fator1),e=length(levels(Fator2)))
                    f2=rep(unique(as.character(Fator2)),length(levels(Fator2)))
                    f1=factor(f1,levels = unique(f1))
                    f2=factor(f2,levels = unique(f2))
                    media=tapply(resp,paste(Fator1,Fator2), mean, na.rm=TRUE)[unique(paste(f1,f2))]
                    desvio=tapply(resp,paste(Fator1,Fator2), sd, na.rm=TRUE)[unique(paste(f1,f2))]
                    graph=data.frame(f1=f1,
                                     f2=f2,
                                     media,
                                     desvio,
                                     letra,letra1,
                                     numero=format(media,digits = dec))
                    numero=paste(graph$numero,graph$letra,graph$letra1,sep="")
                    graph$numero=numero
                    colint=ggplot(graph, aes(x=f2,
                                             y=media,
                                             fill=f1))+
                        ylab(ylab)+xlab(xlab)+
                        theme+
                        labs(fill=fac.names[1])+
                        geom_col(position = "dodge",color="black")+
                        geom_errorbar(aes(ymin=media-desvio,
                                          ymax=media+desvio),
                                      width=0.3,position = position_dodge(width=0.9))+
                        geom_text(aes(y=media+sup+if(sup<0){-desvio}else{desvio},
                                      label=numero),
                                  position = position_dodge(width=0.9),angle=angle.label, hjust=hjust)+
                        theme(text=element_text(size=12),
                              axis.text = element_text(size=12,color="black"),
                              axis.title = element_text(size=12,color="black"))
                    colint1=colint
                    print(colint)
                    letras=paste(graph$letra,graph$letra1,sep="")
                    matriz=data.frame(t(matrix(paste(format(graph$media,digits = dec),letras),ncol = length(levels(Fator1)))))
                    rownames(matriz)=levels(Fator1)
                    colnames(matriz)=levels(Fator2)
                    cat(green(bold("\n------------------------------------------\n")))
                    cat(green(bold("Final table")))
                    cat(green(bold("\n------------------------------------------\n")))
                    print(matriz)
                    cat("\n\nAverages followed by the same lowercase letter in the column and \nuppercase in the row do not differ by the",mcomp,"(p<",alpha.t,")")
                }
        if(quali[1]==FALSE | quali[2]==FALSE){
            if(quali[1]==FALSE){
                if (mcomp == "tukey"){
                    for (i in 1:nv1) {
                        trati=as.factor(fatores[, 2][Fator1 == lf1[i]])
                        respi=resp[Fator1 == lf1[i]]
                        mod=aov(respi~trati)
                        tukey=TUKEY(mod,"trati",anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                        if(transf !="1"){tukey$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(tukey$groups)]}
                        cat("\n----------------------\n")
                        cat("Multiple comparison of F2 within level",lf1[i],"of F1")
                        cat("\n----------------------\n")
                        print(tukey$groups)}}
                if (mcomp == "duncan"){
                    for (i in 1:nv1) {
                        trati=as.factor(fatores[, 2][Fator1 == lf1[i]])
                        respi=resp[Fator1 == lf1[i]]
                        duncan=duncan(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                        if(transf !="1"){duncan$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(duncan$groups)]}
                        cat("\n----------------------\n")
                        cat("Multiple comparison of F2 within level",lf1[i],"of F1")
                        cat("\n----------------------\n")
                        print(duncan$groups)}}
                if (mcomp == "lsd"){
                    for (i in 1:nv1) {
                        trati=as.factor(fatores[, 2][Fator1 == lf1[i]])
                        respi=resp[Fator1 == lf1[i]]
                        lsd=LSD(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                        if(transf !="1"){lsd$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(lsd$groups)]}
                        cat("\n----------------------\n")
                        cat("Multiple comparison of F2 within level",lf1[i],"of F1")
                        cat("\n----------------------\n")
                        print(lsd$groups)}}
                if (mcomp == "sk"){
                    for (i in 1:nv1) {
                        trati=fatores[, 2][Fator1 == lf1[i]]
                        trati=factor(trati,levels = unique(trati))
                        respi=resp[Fator1 == lf1[i]]
                        sk=sk(respi,trati,anavaF3$Df[9],anavaF3$`Sum Sq`[9],alpha.t)
                        if(transf !="1"){sk$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(sk)]}
                        cat("\n----------------------\n")
                        cat("Multiple comparison of F2 within level",lf1[i],"of F1")
                        cat("\n----------------------\n")
                        print(sk)}}}
            if(quali[1]==FALSE){
                Fator1a=fator1a#as.numeric(as.character(Fator1))
                colint1=polynomial2(Fator1a,
                                    response,
                                    Fator3,
                                    grau = grau,
                                    ylab=ylab,
                                    xlab=xlab,
                                    theme=theme)}
            if(quali[2]==FALSE){
                if (mcomp == "tukey"){
                    for (i in 1:nv2) {
                        trati=fatores[, 1][Fator2 == lf2[i]]
                        trati=factor(trati,levels = unique(trati))
                        respi=resp[Fator2 == lf2[i]]
                        mod=aov(respi~trati)
                        tukey=TUKEY(mod,"trati",anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                        cat("\n----------------------\n")
                        cat("Multiple comparison of F1 within level",lf2[i],"of F2")
                        cat("\n----------------------\n")
                        print(tukey$groups)}}
                if (mcomp == "duncan"){
                    for (i in 1:nv2) {
                        trati=fatores[, 1][Fator2 == lf2[i]]
                        trati=factor(trati,levels = unique(trati))
                        respi=resp[Fator2 == lf2[i]]
                        duncan=duncan(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                        if(transf !="1"){duncan$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(duncan$groups)]}
                        cat("\n----------------------\n")
                        cat("Multiple comparison of F1 within level",lf2[i],"of F2")
                        cat("\n----------------------\n")
                        print(duncan$groups)}}
                if (mcomp == "lsd"){
                    for (i in 1:nv2) {
                        trati=fatores[, 1][Fator2 == lf2[i]]
                        trati=factor(trati,levels = unique(trati))
                        respi=resp[Fator2 == lf2[i]]
                        lsd=LSD(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                        if(transf !="1"){lsd$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(lsd$groups)]}
                        cat("\n----------------------\n")
                        cat("Multiple comparison of F1 within level",lf2[i],"of F2")
                        cat("\n----------------------\n")
                        print(lsd$groups)}}
                if (mcomp == "sk"){
                    for (i in 1:nv2) {
                        trati=fatores[, 1][Fator2 == lf2[i]]
                        trati=factor(trati,levels = unique(trati))
                        respi=resp[Fator2 == lf2[i]]
                        sk=sk(respi,trati,anavaF3$Df[9],anavaF3$`Sum Sq`[9],alpha.t)
                        if(transf !="1"){sk$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(sk$groups)]}
                        cat("\n----------------------\n")
                        cat("Multiple comparison of F1 within level",lf2[i],"of F2")
                        cat("\n----------------------\n")
                        print(sk)}}
                }
            if(quali[2]==FALSE){
                Fator2a=fator2a#as.numeric(as.character(Fator2))
                colint1=polynomial2(Fator2a,
                                    response,
                                    Fator1,
                                    grau = grau,
                                    ylab=ylab,
                                    xlab=xlab,
                                    theme=theme)}
            cat(green("To edit graphical parameters, I suggest analyzing using the \"polynomial2\" command\n"))}

        #Checar o Fator3
        if(anavaF3[6,5]>alpha.f && anavaF3[7,5]>alpha.f) {


                i<-3
                {
                    #Para os fatores QUALITATIVOS, teste de Tukey
                    if(quali[i]==TRUE && anavaF3[i,5]<=alpha.f) {
                        cat(green(bold("\n------------------------------------------\n")))
                        cat(green(italic('Analyzing the simple effects of the factor ',fac.names[i])))
                        cat(green(bold("\n------------------------------------------\n")))
                        cat(fac.names[i])
                        if(mcomp=='tukey'){letra=TUKEY(resp,fatores[,i],anavaF3[9,1],anavaF3[9,3],alpha.t)
                        letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
                        if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
                        if(mcomp=="sk"){
                            ad=data.frame(Fator1,Fator2,Fator3)
                            letra=SK(anava,colnames(ad[i]))
                            letra1=data.frame(resp=letra$m.inf[,1],groups=letters[letra$groups])
                            letra1$resp=as.numeric(as.character(letra1$resp))
                            if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
                        if(mcomp=="duncan"){
                            ad=data.frame(Fator1,Fator2,Fator3)
                            letra <- duncan(anava, colnames(ad[i]), alpha=alpha.t)
                            letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
                            if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
                        if(mcomp=="lsd"){
                            ad=data.frame(Fator1,Fator2,Fator3)
                            letra <- LSD(anava, colnames(ad[i]), alpha=alpha.t)
                            letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
                            if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
                        print(letra1)
                        cat(green(bold("\n------------------------------------------\n")))
                        dadosm=data.frame(letra1,
                                          media=tapply(response, c(fatores[i]), mean, na.rm=TRUE)[rownames(letra1)],
                                          desvio=tapply(response, c(fatores[i]), sd, na.rm=TRUE)[rownames(letra1)])
                        dadosm$Tratamentos=factor(rownames(dadosm),levels = unique(unlist(fatores[i])))
                        dadosm$limite=dadosm$media+dadosm$desvio
                        dadosm=dadosm[as.character(unique(unlist(fatores[i]))),]

                        if(addmean==TRUE){dadosm$letra=paste(format(dadosm$media,digits = dec),dadosm$groups)}
                        if(addmean==FALSE){dadosm$letra=dadosm$groups}
                        media=dadosm$media
                        desvio=dadosm$desvio
                        Tratamentos=dadosm$Tratamentos
                        letra=dadosm$letra

                        grafico=ggplot(dadosm,
                                       aes(x=Tratamentos,
                                           y=media))
                        if(fill=="trat"){grafico=grafico+
                            geom_col(aes(fill=Tratamentos),color=1)}
                        else{grafico=grafico+
                            geom_col(aes(fill=Tratamentos),fill=fill,color=1)}
                        if(errorbar==TRUE){grafico=grafico+
                            geom_text(aes(y=media+sup+if(sup<0){-desvio}else{desvio},
                                          label=letra),family=family,angle=angle.label, hjust=hjust)}
                        if(errorbar==FALSE){grafico=grafico+
                            geom_text(aes(y=media+sup,label=letra),family=family,angle=angle.label, hjust=hjust)}
                        if(errorbar==TRUE){grafico=grafico+
                            geom_errorbar(data=dadosm,
                                          aes(ymin=media-desvio,
                                              ymax=media+desvio,color=1), color="black",width=0.3)}
                        grafico1=grafico+theme+
                            ylab(ylab)+
                            xlab(xlab)+
                            theme(text = element_text(size=textsize,color="black", family = family),
                                  axis.text = element_text(size=textsize,color="black", family = family),
                                  axis.title = element_text(size=textsize,color="black", family = family),
                                  legend.position = "none")
                        print(grafico1)}

                        }

                    #Para os fatores QUANTITATIVOS, regressao
                    if(quali[i]==FALSE && anavaF3[i,5]<=alpha.f){
                        cat(green(bold("\n------------------------------------------\n")))
                        cat('Analyzing the simple effects of the factor ',fac.names[3])
                        cat(green(bold("\n------------------------------------------\n")))
                        cat(fac.names[i])
                        grafico1=polynomial(resp, fatores[,i])
                        cat(green("To edit graphical parameters, I suggest analyzing using the \"polynomial\" command"))}
                }
            }

    #####################################################################################################
    #Interacao Fator1*Fator3       + fator2
    #####################################################################################################
    if(anavaF3[8,5]>alpha.f && anavaF3[6,5]<=alpha.f){
    cat(green(bold("\n------------------------------------------\n")))
        cat(green(bold("Interaction",paste(fac.names[1],'*',fac.names[3],sep='')," significant: unfolding the interaction")))
    cat(green(bold("\n------------------------------------------")))

        #Desdobramento de FATOR 1 dentro do niveis de FATOR 3
        cat(green(bold("\n------------------------------------------\n")))
        cat("Analyzing ", fac.names[1], ' within the combination of levels ', fac.names[3])
        cat(green(bold("\n------------------------------------------\n")))
        des<-aov(resp~Fator3/Fator1+Fator2+Fator3+Fator2:Fator3+Fator1:Fator2:Fator3)
        l<-vector('list',nv3)
        names(l)<-names(summary(Fator3))
        v<-numeric(0)
        for(j in 1:nv3) {
            for(i in 0:(nv1-2)) v<-cbind(v,i*nv3+j)
            l[[j]]<-v
            v<-numeric(0)
        }
        des1<-summary(des,split=list('Fator3:Fator1'=l))[[1]]
        des1a=des1[-c(1,2,3,length(des1[,1]),length(des1[,1])-1,length(des1[,1])-2),]
        print(des1a)

        # Teste de Tukey
        if(quali[1]==TRUE & quali[3]==TRUE){
        if (mcomp == "tukey"){
            tukeygrafico=c()
            ordem=c()
            for (i in 1:nv3) {
                trati=fatores[, 1][Fator3 == lf3[i]]
                trati=factor(trati,levels = unique(trati))
                respi=resp[Fator3 == lf3[i]]
                tukey=TUKEY(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                tukeygrafico[[i]]=tukey$groups[levels(trati),2]
                ordem[[i]]=rownames(tukey$groups[levels(trati),])
                }
            letra=unlist(tukeygrafico)
            datag=data.frame(letra,ordem=unlist(ordem))
            datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
            datag=datag[order(datag$ordem),]
            letra=datag$letra}
        if (mcomp == "duncan"){
            duncangrafico=c()
            ordem=c()
            for (i in 1:nv3) {
                trati=fatores[, 1][Fator3 == lf3[i]]
                trati=factor(trati,levels = unique(trati))
                respi=resp[Fator3 == lf3[i]]
                duncan=duncan(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                duncangrafico[[i]]=duncan$groups[levels(trati),2]
                ordem[[i]]=rownames(duncan$groups[levels(trati),])
                }
            letra=unlist(duncangrafico)
            datag=data.frame(letra,ordem=unlist(ordem))
            datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
            datag=datag[order(datag$ordem),]
            letra=datag$letra}
        if (mcomp == "lsd"){
            lsdgrafico=c()
            ordem=c()
            for (i in 1:nv3) {
                trati=fatores[, 1][Fator3 == lf3[i]]
                trati=factor(trati,levels = unique(trati))
                respi=resp[Fator3 == lf3[i]]
                lsd=LSD(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                lsdgrafico[[i]]=lsd$groups[levels(trati),2]
                ordem[[i]]=rownames(lsd$groups[levels(trati),])
                }
            letra=unlist(lsdgrafico)
            datag=data.frame(letra,ordem=unlist(ordem))
            datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
            datag=datag[order(datag$ordem),]
            letra=datag$letra}
        if (mcomp == "sk"){
            skgrafico=c()
            ordem=c()
            for (i in 1:nv3) {
                trati=fatores[, 1][Fator3 == lf3[i]]
                trati=factor(trati,levels = unique(trati))
                respi=resp[Fator3 == lf3[i]]
                sk=sk(respi,trati,anavaF3$Df[9],anavaF3$`Sum Sq`[9],alpha.t)
                skgrafico[[i]]=sk[levels(trati),2]
                ordem[[i]]=rownames(sk[levels(trati),])
                }
            letra=unlist(skgrafico)
            datag=data.frame(letra,ordem=unlist(ordem))
            datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
            datag=datag[order(datag$ordem),]
            letra=datag$letra}}

        # Desdobramento de F3 dentro de F1

        cat(green(bold("\n------------------------------------------\n")))
        cat("Analyzing ", fac.names[3], " inside of the level of ",fac.names[1])
        cat(green(bold("\n------------------------------------------\n")))
        des<-aov(resp~Fator1/Fator3+Fator1+Fator2+Fator2:Fator1+Fator1:Fator2:Fator3)
        l<-vector('list',nv1)
        names(l)<-names(summary(Fator1))
        v<-numeric(0)
        for(j in 1:nv1) {
                for(i in 0:(nv3-2)) v<-cbind(v,i*nv1+j)
                l[[j]]<-v
                v<-numeric(0)
            }
        des1<-summary(des,split=list('Fator1:Fator3'=l))[[1]]
        des1a=des1[-c(1,2,3,length(des1[,1]),length(des1[,1])-1,length(des1[,1])-2),]
        print(des1a)

        #-------------------------------------
        # Teste de Tukey
        #-------------------------------------
        if(quali[1]==TRUE & quali[3]==TRUE){
        if (mcomp == "tukey"){
                    tukeygrafico1=c()
                    for (i in 1:nv1) {
                        trati=fatores[, 3][Fator1 == lf1[i]]
                        trati=factor(trati,levels = unique(trati))
                        respi=resp[Fator1 == lf1[i]]
                        mod=aov(respi~trati)
                        tukey=TUKEY(mod,"trati",anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                        tukeygrafico1[[i]]=tukey$groups[levels(trati),2]
                        }
                    letra1=unlist(tukeygrafico1)
                    letra1=toupper(letra1)}
        if (mcomp == "duncan"){
                    duncangrafico=c()
                    ordem=c()
                    for (i in 1:nv3) {
                        trati=fatores[, 3][Fator1 == lf1[i]]
                        trati=factor(trati,levels = unique(trati))
                        respi=resp[Fator1 == lf1[i]]
                        duncan=duncan(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                        duncangrafico[[i]]=duncan$groups[levels(trati),2]
                        ordem[[i]]=rownames(duncan$groups[levels(trati),])
                        }
                    letra=unlist(duncangrafico)
                    datag=data.frame(letra,ordem=unlist(ordem))
                    datag=datag[order(datag$ordem),]
                    letra=datag$letra}
        if (mcomp == "lsd"){
                    lsdgrafico=c()
                    ordem=c()
                    for (i in 1:nv3) {
                        trati=fatores[, 3][Fator1 == lf1[i]]
                        trati=factor(trati,levels = unique(trati))
                        respi=resp[Fator1 == lf1[i]]
                        lsd=LSD(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                        lsdgrafico[[i]]=lsd$groups[levels(trati),2]
                        ordem[[i]]=rownames(lsd$groups[levels(trati),])
                        }
                    letra=unlist(lsdgrafico)
                    datag=data.frame(letra,ordem=unlist(ordem))
                    datag=datag[order(datag$ordem),]
                    letra=datag$letra}
        if (mcomp == "sk"){
                    skgrafico1=c()
                    for (i in 1:nv1) {
                        trati=fatores[, 3][Fator1 == lf1[i]]
                        trati=factor(trati,levels = unique(trati))
                        respi=resp[Fator1 == lf1[i]]
                        mod=aov(respi~trati)
                        sk=sk(respi,trati,anavaF3$Df[9],anavaF3$`Sum Sq`[9],alpha.t)
                        skgrafico1[[i]]=sk[levels(trati),2]
                        }
                    letra1=unlist(skgrafico1)
                    letra1=toupper(letra1)}}

        # -----------------------------
        # Gráfico de colunas
        #------------------------------
        if(quali[1] & quali[3]==TRUE){
        f1=rep(levels(Fator1),e=length(levels(Fator3)))
        f3=rep(unique(as.character(Fator3)),length(levels(Fator1)))
        f1=factor(f1,levels = unique(f1))
        f3=factor(f3,levels = unique(f3))
        media=tapply(response,paste(Fator1,Fator3), mean, na.rm=TRUE)[unique(paste(f1,f3))]
        desvio=tapply(response,paste(Fator1,Fator3), sd, na.rm=TRUE)[unique(paste(f1,f3))]
        graph=data.frame(f1=f1,
                         f3=f3,
                         media,
                         desvio,
                         letra,
                         letra1,
                         numero=format(media,digits = dec))
        numero=paste(graph$numero,graph$letra,graph$letra1,sep="")
        graph$numero=numero
        colint=ggplot(graph,
                      aes(x=f3,
                          y=media,
                          fill=f1))+
            geom_col(position = "dodge",color="black")+
            ylab(ylab)+
            xlab(xlab)+
            theme+
            labs(fill=fac.names[1])+
            geom_errorbar(aes(ymin=media-desvio,
                              ymax=media+desvio),
                          width=0.3,
                          position = position_dodge(width=0.9))+
            geom_text(aes(y=media+desvio+sup,
                          label=numero),
                      position = position_dodge(width=0.9),angle=angle.label, hjust=hjust)+
            theme(text=element_text(size=12),
                  axis.text = element_text(size=12,color="black"),
                  axis.title = element_text(size=12,color="black"))
        colint2=colint
        print(colint)
        letras=paste(graph$letra,graph$letra1,sep="")
        matriz=data.frame(t(matrix(paste(format(graph$media,digits = dec),letras),ncol = length(levels(Fator1)))))
        rownames(matriz)=levels(Fator1)
        colnames(matriz)=levels(Fator2)
        cat(green(bold("\n------------------------------------------\n")))
        cat(green(bold("Final table")))
        cat(green(bold("\n------------------------------------------\n")))
        print(matriz)
        cat("\n\nAverages followed by the same lowercase letter in the column and \nuppercase in the row do not differ by the",mcomp,"(p<",alpha.t,")")
        }
        if(quali[1]==FALSE | quali[3]==FALSE){
            if(quali[1]==FALSE){
                if (mcomp == "tukey"){
                    for (i in 1:nv1) {
                        trati=fatores[, 3][Fator1 == lf1[i]]
                        trati=factor(trati,levels = unique(trati))
                        respi=resp[Fator1 == lf1[i]]
                        mod=aov(respi~trati)
                        tukey=TUKEY(mod,"trati",anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                        if(transf !="1"){tukey$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(tukey$groups)]}
                        cat("\n----------------------\n")
                        cat("Multiple comparison of F3 within level",lf1[i],"of F1")
                        cat("\n----------------------\n")
                        print(tukey$groups)}}
                if (mcomp == "duncan"){
                    for (i in 1:nv3) {
                        trati=fatores[, 3][Fator1 == lf1[i]]
                        trati=factor(trati,levels = unique(trati))
                        respi=resp[Fator1 == lf1[i]]
                        duncan=duncan(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                        if(transf !="1"){duncan$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(duncan$groups)]}
                        cat("\n----------------------\n")
                        cat("Multiple comparison of F3 within level",lf1[i],"of F1")
                        cat("\n----------------------\n")
                        print(duncan$groups)}}
                if (mcomp == "lsd"){
                    for (i in 1:nv3) {
                        trati=fatores[, 3][Fator1 == lf1[i]]
                        trati=factor(trati,levels = unique(trati))
                        respi=resp[Fator1 == lf1[i]]
                        lsd=LSD(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                        if(transf !="1"){lsd$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(lsd$groups)]}
                        cat("\n----------------------\n")
                        cat("Multiple comparison of F3 within level",lf1[i],"of F1")
                        cat("\n----------------------\n")
                        print(lsd$groups)}}
                if (mcomp == "sk"){
                    for (i in 1:nv1) {
                        trati=fatores[, 3][Fator1 == lf1[i]]
                        trati=factor(trati,levels = unique(trati))
                        respi=resp[Fator1 == lf1[i]]
                        mod=aov(respi~trati)
                        sk=sk(respi,trati,anavaF3$Df[8],anavaF3$`Mean Sq`[8],alpha.t)
                        if(transf !="1"){sk$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(sk)]}
                        cat("\n----------------------\n")
                        cat("Multiple comparison of F3 within level",lf1[i],"of F1")
                        cat("\n----------------------\n")
                        print(sk)}}}
            if(quali[1]==FALSE){
                Fator1a=fator1a#as.numeric(as.character(Fator3))
                colint2=polynomial2(Fator1a,
                                    response,
                                    Fator3,
                                    grau = grau,
                                    ylab=ylab,
                                    xlab=xlab,
                                    theme=theme)}
            if(quali[3]==FALSE){
                if (mcomp == "tukey"){
                    for (i in 1:nv3) {
                        trati=fatores[, 1][Fator3 == lf3[i]]
                        trati=factor(trati,levels = unique(trati))
                        respi=resp[Fator3 == lf3[i]]
                        tukey=TUKEY(respi,trati,anavaF3$Df[9],anavaF3$`Sum Sq`[9],alpha.t)
                        if(transf !="1"){tukey$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(tukey$groups)]}
                        cat("\n----------------------\n")
                        cat("Multiple comparison of F1 within level",lf3[i],"of F3")
                        cat("\n----------------------\n")
                        print(tukey$groups)}}
                if (mcomp == "duncan"){
                    for (i in 1:nv3) {
                        trati=fatores[, 1][Fator3 == lf3[i]]
                        trati=factor(trati,levels = unique(trati))
                        respi=resp[Fator3 == lf3[i]]
                        duncan=duncan(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                        if(transf !="1"){duncan$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(duncan$groups)]}
                        cat("\n----------------------\n")
                        cat("Multiple comparison of F1 within level",lf3[i],"of F3")
                        cat("\n----------------------\n")
                        print(duncan$groups)}}
                if (mcomp == "lsd"){
                    for (i in 1:nv3) {
                        trati=fatores[, 1][Fator3 == lf3[i]]
                        trati=factor(trati,levels = unique(trati))
                        respi=resp[Fator3 == lf3[i]]
                        lsd=LSD(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                        if(transf !="1"){lsd$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(lsd$groups)]}
                        cat("\n----------------------\n")
                        cat("Multiple comparison of F1 within level",lf3[i],"of F3")
                        cat("\n----------------------\n")
                        print(lsd$groups)}}
                if (mcomp == "sk"){
                    for (i in 1:nv3) {
                        trati=fatores[, 1][Fator3 == lf3[i]]
                        trati=factor(trati,levels = unique(trati))
                        respi=resp[Fator3 == lf3[i]]
                        sk=sk(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                        if(transf !="1"){sk$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(sk)]}
                        cat("\n----------------------\n")
                        cat("Multiple comparison of F1 within level",lf3[i],"of F3")
                        cat("\n----------------------\n")
                        print(sk)}}}
            if(quali[3]==FALSE){
                Fator3a=fator3a#as.numeric(as.character(Fator3))
                colint2=polynomial2(Fator3a,
                                    response,
                                    Fator1,
                                    grau = grau,
                                    ylab=ylab,
                                    xlab=xlab,
                                    theme=theme)}
        cat(green("To edit graphical parameters, I suggest analyzing using the \"polynomial2\" command\n"))
        }

        #Checar o Fator2
        if(anavaF3[5,5]>alpha.f && anavaF3[7,5]>alpha.f) {


                i<-2
                {
                    #Para os fatores QUALITATIVOS, teste de Tukey
                    if(quali[i]==TRUE && anavaF3[i,5]<=alpha.f) {
                        cat(green(bold("\n------------------------------------------\n")))
                        cat(green(italic('Analyzing the simple effects of the factor ',fac.names[i])))
                        cat(green(bold("\n------------------------------------------\n")))
                        cat(fac.names[i])
                        if(mcomp=='tukey'){letra=TUKEY(resp,fatores[,i],anavaF3[9,1],anavaF3[9,3],alpha.t)
                        letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
                        if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
                        if(mcomp=="sk"){
                            ad=data.frame(Fator1,Fator2,Fator3)
                            letra=SK(anava,colnames(ad[i]))
                            letra1=data.frame(resp=letra$m.inf[,1],groups=letters[letra$groups])
                            letra1$resp=as.numeric(as.character(letra1$resp))
                            if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
                        if(mcomp=="duncan"){
                            ad=data.frame(Fator1,Fator2,Fator3)
                            letra <- duncan(anava, colnames(ad[i]), alpha=alpha.t)
                            letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
                            if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
                        if(mcomp=="lsd"){
                            ad=data.frame(Fator1,Fator2,Fator3)
                            letra <- LSD(anava, colnames(ad[i]), alpha=alpha.t)
                            letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
                            if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
                        print(letra1)
                        cat(green(bold("\n------------------------------------------")))
                        dadosm=data.frame(letra1,
                                          media=tapply(response, c(fatores[i]), mean, na.rm=TRUE)[rownames(letra1)],
                                          desvio=tapply(response, c(fatores[i]), sd, na.rm=TRUE)[rownames(letra1)])
                        dadosm$Tratamentos=factor(rownames(dadosm),levels = unique(unlist(fatores[i])))
                        dadosm$limite=dadosm$media+dadosm$desvio
                        dadosm=dadosm[as.character(unique(unlist(fatores[i]))),]
                        if(addmean==TRUE){dadosm$letra=paste(format(dadosm$media,digits = dec),dadosm$groups)}
                        if(addmean==FALSE){dadosm$letra=dadosm$groups}
                        media=dadosm$media
                        desvio=dadosm$desvio
                        Tratamentos=dadosm$Tratamentos
                        letra=dadosm$letra

                        grafico=ggplot(dadosm,
                                       aes(x=Tratamentos,
                                           y=media))
                        if(fill=="trat"){grafico=grafico+
                            geom_col(aes(fill=Tratamentos),color=1)}
                        else{grafico=grafico+
                            geom_col(aes(fill=Tratamentos),fill=fill,color=1)}
                        if(errorbar==TRUE){grafico=grafico+
                            geom_text(aes(y=media+sup+if(sup<0){-desvio}else{desvio},label=letra),family=family,angle=angle.label, hjust=hjust)}
                        if(errorbar==FALSE){grafico=grafico+
                            geom_text(aes(y=media+sup,label=letra),family=family,angle=angle.label, hjust=hjust)}
                        if(errorbar==TRUE){grafico=grafico+
                            geom_errorbar(data=dadosm,aes(ymin=media-desvio,
                                                          ymax=media+desvio,color=1),
                                          color="black",width=0.3)
                        grafico2=grafico+theme+
                            ylab(ylab)+
                            xlab(xlab)+
                            theme(text = element_text(size=textsize,color="black", family = family),
                                  axis.text = element_text(size=textsize,color="black", family = family),
                                  axis.title = element_text(size=textsize,color="black", family = family),
                                  legend.position = "none")
                        print(grafico2)}
                    }

                    #Para os fatores QUANTITATIVOS, regressao
                    if(quali[i]==FALSE && anavaF3[i,5]<=alpha.f){
                        cat(green(bold("\n------------------------------------------\n")))
                        cat('Analyzing the simple effects of the factor ',fac.names[2])
                        cat(green(bold("\n------------------------------------------\n")))
                        cat(fac.names[i])
                        grafico2=polynomial(resp, fatores[,i])
                        cat(green("To edit graphical parameters, I suggest analyzing using the \"polynomial\" command"))
                    }

                    cat('\n')
                }
            }
    }

    ######################################################################################################################
    #Interacao Fator2*Fator3     + fator1
    ######################################################################################################################
    if(anavaF3[8,5]>alpha.f && anavaF3[7,5]<=alpha.f){
        cat(green(bold("\n------------------------------------------\n")))
        cat(green(bold("Interaction",paste(fac.names[2],'*',fac.names[3],sep='')," significant: unfolding the interaction")))
        cat(green(bold("\n------------------------------------------\n")))
        #Desdobramento de FATOR 2 dentro do niveis de FATOR 3
        cat(green(bold("\n------------------------------------------\n")))
        cat("Analyzing ", fac.names[2], ' within the combination of levels ', fac.names[3])
        cat("\n-------------------------------------------------\n")
        des<-aov(resp~Fator3/Fator2+Fator1+Fator3+Fator1:Fator3+Fator1:Fator2:Fator3)
        l<-vector('list',nv3)
        names(l)<-names(summary(Fator3))
        v<-numeric(0)
        for(j in 1:nv3) {
            for(i in 0:(nv2-2)) v<-cbind(v,i*nv3+j)
            l[[j]]<-v
            v<-numeric(0)
        }
        des1<-summary(des,split=list('Fator3:Fator2'=l))[[1]]
        des1a=des1[-c(1,2,3,length(des1[,1]),length(des1[,1])-1,length(des1[,1])-2),]
        print(des1a)

        # Teste de Tukey
        if(quali[2]==TRUE & quali[3]==TRUE){
            if (mcomp == "tukey"){
                tukeygrafico=c()
                ordem=c()
                for (i in 1:nv3) {
                    trati=fatores[, 2][Fator3 == lf3[i]]
                    trati=factor(trati,levels = unique(trati))
                    respi=resp[Fator3 == lf3[i]]
                    mod=aov(respi~trati)
                    tukey=TUKEY(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                    tukeygrafico[[i]]=tukey$groups[levels(trati),2]
                    ordem[[i]]=rownames(tukey$groups[levels(trati),])
                    }
                letra=unlist(tukeygrafico)
                datag=data.frame(letra,ordem=unlist(ordem))
                datag=data.frame(letra,ordem=unlist(ordem))
                datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
                datag=datag[order(datag$ordem),]
                letra=datag$letra}
            if (mcomp == "duncan"){
                duncangrafico=c()
                ordem=c()
                for (i in 1:nv3) {
                    trati=fatores[, 2][Fator3 == lf3[i]]
                    trati=factor(trati,levels = unique(trati))
                    respi=resp[Fator3 == lf3[i]]
                    duncan=duncan(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                    duncangrafico[[i]]=duncan$groups[levels(trati),2]
                    ordem[[i]]=rownames(duncan$groups[levels(trati),])
                    }
                letra=unlist(duncangrafico)
                datag=data.frame(letra,ordem=unlist(ordem))
                datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
                datag=datag[order(datag$ordem),]
                letra=datag$letra}
            if (mcomp == "lsd"){
                lsdgrafico=c()
                ordem=c()
                for (i in 1:nv3) {
                    trati=fatores[, 2][Fator3 == lf3[i]]
                    trati=factor(trati,levels = unique(trati))
                    respi=resp[Fator3 == lf3[i]]
                    lsd=LSD(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                    lsdgrafico[[i]]=lsd$groups[levels(trati),2]
                    ordem[[i]]=rownames(lsd$groups[levels(trati),])
                    }
                letra=unlist(lsdgrafico)
                datag=data.frame(letra,ordem=unlist(ordem))
                datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
                datag=datag[order(datag$ordem),]
                letra=datag$letra}
            if (mcomp == "sk"){
                skgrafico=c()
                ordem=c()
                for (i in 1:nv3) {
                    trati=fatores[, 2][Fator3 == lf3[i]]
                    trati=factor(trati,levels = unique(trati))
                    respi=resp[Fator3 == lf3[i]]
                    sk=sk(respi,trati,anavaF3$Df[9],anavaF3$`Sum Sq`[9],alpha.t)
                    skgrafico[[i]]=sk[levels(trati),2]
                    ordem[[i]]=rownames(sk[levels(trati),])
                    }
                letra=unlist(skgrafico)
                datag=data.frame(letra,ordem=unlist(ordem))
                datag$ordem=factor(datag$ordem,levels = unique(datag$ordem))
                datag=datag[order(datag$ordem),]
                letra=datag$letra}}

        # Desdobramento de F2 dentro de F1

            cat(green(bold("\n------------------------------------------\n")))
            cat("Analyzing ", fac.names[3], " inside of the level of ",fac.names[2])
            cat(green(bold("\n------------------------------------------\n")))
            cat("\n")
            des<-aov(resp~Fator2/Fator3+Fator1+Fator2+Fator1:Fator2+Fator1:Fator2:Fator3)
            l<-vector('list',nv2)
            names(l)<-names(summary(Fator2))
            v<-numeric(0)
            for(j in 1:nv2) {
                for(i in 0:(nv3-2)) v<-cbind(v,i*nv2+j)
                l[[j]]<-v
                v<-numeric(0)
            }
            des1<-summary(des,split=list('Fator2:Fator3'=l))[[1]]
            des1a=des1[-c(1,2,3,length(des1[,1]),length(des1[,1])-1,length(des1[,1])-2),]
            print(des1a)

            #-------------------------------------
            # Teste de Tukey
            #-------------------------------------
            if(quali[2]==TRUE & quali[3]==TRUE){
                if (mcomp == "tukey"){
                    tukeygrafico1=c()
                    for (i in 1:nv2) {
                        trati=fatores[, 3][Fator2 == lf2[i]]
                        trati=factor(trati,levels = unique(trati))
                        respi=resp[Fator2 == lf2[i]]
                        mod=aov(respi~trati)
                        tukey=TUKEY(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                        tukeygrafico1[[i]]=tukey$groups[levels(trati),2]
                        }
                    letra1=unlist(tukeygrafico1)
                    letra1=toupper(letra1)}
                if (mcomp == "duncan"){
                    duncangrafico1=c()
                    for (i in 1:nv2) {
                        trati=fatores[, 3][Fator2 == lf2[i]]
                        trati=factor(trati,levels = unique(trati))
                        respi=resp[Fator2 == lf2[i]]
                        duncan=duncan(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                        duncangrafico1[[i]]=duncan$groups[levels(trati),2]
                        }
                    letra1=unlist(duncangrafico1)
                    letra1=toupper(letra1)}
                if (mcomp == "lsd"){
                    lsdgrafico1=c()
                    for (i in 1:nv2) {
                        trati=fatores[, 3][Fator2 == lf2[i]]
                        trati=factor(trati,levels = unique(trati))
                        respi=resp[Fator2 == lf2[i]]
                        lsd=LSD(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                        lsdgrafico1[[i]]=lsd$groups[levels(trati),2]
                        }
                    letra1=unlist(lsdgrafico1)
                    letra1=toupper(letra1)}
                if (mcomp == "sk"){
                    skgrafico1=c()
                    for (i in 1:nv2) {
                        trati=fatores[, 3][Fator2 == lf2[i]]
                        trati=factor(trati,levels = unique(trati))
                        respi=resp[Fator2 == lf2[i]]
                        sk=sk(respi,trati,anavaF3$Df[9],anavaF3$`Sum Sq`[9],alpha.t)
                        skgrafico1[[i]]=sk[levels(trati),2]
                        }
                    letra1=unlist(skgrafico1)
                    letra1=toupper(letra1)}}

            # -----------------------------
            # Gráfico de colunas
            #------------------------------
            if(quali[2] & quali[3]==TRUE){
                f2=rep(levels(Fator2),e=length(levels(Fator3)))
                f3=rep(unique(as.character(Fator3)),length(levels(Fator2)))
                f2=factor(f2,levels = unique(f2))
                f3=factor(f3,levels = unique(f3))
                media=tapply(response,paste(Fator2,Fator3), mean, na.rm=TRUE)[unique(paste(f2,f3))]
                desvio=tapply(response,paste(Fator2,Fator3), sd, na.rm=TRUE)[unique(paste(f2,f3))]
                graph=data.frame(f2=f2,
                                 f3=f3,
                                 media,
                                 desvio,
                                 letra,letra1,
                                 numero=format(media,digits = dec))
                numero=paste(graph$numero,graph$letra,graph$letra1,sep="")
                graph$numero=numero
                colint=ggplot(graph, aes(x=f3,
                                         y=media,
                                         fill=f2))+
                        geom_col(position = "dodge",color="black")+
                        ylab(ylab)+xlab(xlab)+
                        theme+
                        labs(fill=fac.names[2])+
                        geom_errorbar(aes(ymin=media-desvio,
                                          ymax=media+desvio),
                                      width=0.3,position = position_dodge(width=0.9))+
                        geom_text(aes(y=media+desvio+sup,
                                      label=numero),
                                  position = position_dodge(width=0.9),angle=angle.label, hjust=hjust)+
                        theme(text=element_text(size=12),
                              axis.text = element_text(size=12,color="black"),
                              axis.title = element_text(size=12,color="black"))
                    colint3=colint
                    print(colint)
                    letras=paste(graph$letra,graph$letra1,sep="")
                    matriz=data.frame(t(matrix(paste(format(graph$media,digits = dec),letras),ncol = length(levels(Fator1)))))
                    rownames(matriz)=levels(Fator1)
                    colnames(matriz)=levels(Fator2)
                    cat(green(bold("\n------------------------------------------\n")))
                    cat(green(bold("Final table")))
                    cat(green(bold("\n------------------------------------------\n")))
                    print(matriz)
                    cat("\n\nAverages followed by the same lowercase letter in the column and \nuppercase in the row do not differ by the",mcomp,"(p<",alpha.t,")")
            }
            if(quali[2]==FALSE | quali[3]==FALSE){
                if(quali[2]==FALSE){
                    if (mcomp == "tukey"){
                        for (i in 1:nv2) {
                            trati=fatores[, 3][Fator2 == lf2[i]]
                            trati=factor(trati,levels = unique(trati))
                            respi=resp[Fator2 == lf2[i]]
                            mod=aov(respi~trati)
                            tukey=TUKEY(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                            if(transf !="1"){tukey$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(tukey$groups)]}
                            cat("\n----------------------\n")
                            cat("Multiple comparison of F3 within level",lf2[i],"of F2")
                            cat("\n----------------------\n")
                            print(tukey$groups)}}
                    if (mcomp == "duncan"){
                        for (i in 1:nv2) {
                            trati=fatores[, 3][Fator2 == lf2[i]]
                            trati=factor(trati,levels = unique(trati))
                            respi=resp[Fator2 == lf2[i]]
                            duncan=duncan(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                            if(transf !="1"){duncan$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(duncan$groups)]}
                            cat("\n----------------------\n")
                            cat("Multiple comparison of F3 within level",lf2[i],"of F2")
                            cat("\n----------------------\n")
                            print(duncan$groups)}}
                    if (mcomp == "lsd"){
                        for (i in 1:nv2) {
                            trati=fatores[, 3][Fator2 == lf2[i]]
                            trati=factor(trati,levels = unique(trati))
                            respi=resp[Fator2 == lf2[i]]
                            lsd=LSD(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                            if(transf !="1"){lsd$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(lsd$groups)]}
                            cat("\n----------------------\n")
                            cat("Multiple comparison of F3 within level",lf2[i],"of F2")
                            cat("\n----------------------\n")
                            print(lsd$groups)}}
                    if (mcomp == "sk"){
                        for (i in 1:nv2) {
                            trati=fatores[, 3][Fator2 == lf2[i]]
                            trati=factor(trati,levels = unique(trati))
                            respi=resp[Fator2 == lf2[i]]
                            sk=sk(respi,trati,anavaF3$Df[9],anavaF3$`Sum Sq`[9],alpha.t)
                            if(transf !="1"){sk$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(sk)]}
                            cat("\n----------------------\n")
                            cat("Multiple comparison of F3 within level",lf2[i],"of F2")
                            cat("\n----------------------\n")
                            print(sk)}}}
                if(quali[2]==FALSE){
                    Fator2a=fator2a#as.numeric(as.character(Fator2))
                    colint3=polynomial2(Fator2a,
                                        response,
                                        Fator3,
                                        grau = grau,
                                        ylab=ylab,
                                        xlab=xlab,
                                        theme=theme)}
                if(quali[3]==FALSE){
                    if (mcomp == "tukey"){
                        for (i in 1:nv3) {
                            trati=fatores[, 2][Fator3 == lf3[i]]
                            trati=factor(trati,levels = unique(trati))
                            respi=resp[Fator3 == lf3[i]]
                            mod=aov(respi~trati)
                            tukey=TUKEY(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                            if(transf !="1"){tukey$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(tukey$groups)]}
                            cat("\n----------------------\n")
                            cat("Multiple comparison of F2 within level",lf3[i],"of F3")
                            cat("\n----------------------\n")
                            print(tukey)}}
                    if (mcomp == "duncan"){
                        for (i in 1:nv3) {
                            trati=fatores[, 2][Fator3 == lf3[i]]
                            trati=factor(trati,levels = unique(trati))
                            respi=resp[Fator3 == lf3[i]]
                            duncan=duncan(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                            if(transf !="1"){duncan$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(duncan$groups)]}
                            cat("\n----------------------\n")
                            cat("Multiple comparison of F2 within level",lf3[i],"of F3")
                            cat("\n----------------------\n")
                            print(duncan)}}
                    if (mcomp == "lsd"){
                        for (i in 1:nv3) {
                            trati=fatores[, 2][Fator3 == lf3[i]]
                            trati=factor(trati,levels = unique(trati))
                            respi=resp[Fator3 == lf3[i]]
                            lsd=LSD(respi,trati,anavaF3$Df[9],anavaF3$`Mean Sq`[9],alpha.t)
                            if(transf !="1"){lsd$groups$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(lsd$groups)]}
                            cat("\n----------------------\n")
                            cat("Multiple comparison of F2 within level",lf3[i],"of F3")
                            cat("\n----------------------\n")
                            print(lsd)}}
                    if (mcomp == "sk"){
                        for (i in 1:nv3) {
                            trati=fatores[, 2][Fator3 == lf3[i]]
                            trati=factor(trati,levels = unique(trati))
                            respi=resp[Fator3 == lf3[i]]
                            sk=sk(respi,trati,anavaF3$Df[9],anavaF3$`Sum Sq`[9],alpha.t)
                            if(transf !="1"){sk$respo=tapply(respi,trati,mean, na.rm=TRUE)[rownames(sk)]}
                            cat("\n----------------------\n")
                            cat("Multiple comparison of F2 within level",lf3[i],"of F3")
                            cat("\n----------------------\n")
                            print(sk)}}}
                if(quali[3]==FALSE){
                    Fator3a=fator3a#as.numeric(as.character(Fator3))
                    colint3=polynomial2(Fator3a,
                                        response,
                                        Fator2,
                                        grau = grau,
                                        ylab=ylab,
                                        xlab=xlab,
                                        theme=theme)}

                cat(green("To edit graphical parameters, I suggest analyzing using the \"polynomial2\" command\n"))
            }

            #Checar o Fator1
            if(anavaF3[5,5]>alpha.f && anavaF3[6,5]>alpha.f) {


                i<-1
                {
                    #Para os fatores QUALITATIVOS, teste de Tukey
                    if(quali[i]==TRUE && anavaF3[i,5]<=alpha.f) {
                        cat(green(bold("\n------------------------------------------\n")))
                        cat(green(italic('Analyzing the simple effects of the factor ',fac.names[i])))
                        cat(green(bold("\n------------------------------------------\n")))
                        cat(fac.names[i])
                        if(mcomp=='tukey'){letra=TUKEY(resp,fatores[,i],anavaF3[9,1],anavaF3[9,3],alpha.t)
                        letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
                        if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
                        if(mcomp=="sk"){
                            ad=data.frame(Fator1,Fator2,Fator3)
                            letra=SK(anava,colnames(ad[i]))
                            letra1=data.frame(resp=letra$m.inf[,1],groups=letters[letra$groups])
                            letra1$resp=as.numeric(as.character(letra1$resp))
                            if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
                        if(mcomp=="duncan"){
                            ad=data.frame(Fator1,Fator2,Fator3)
                            letra <- duncan(anava, colnames(ad[i]), alpha=alpha.t)
                            letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
                            if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
                        if(mcomp=="lsd"){
                            ad=data.frame(Fator1,Fator2,Fator3)
                            letra <- LSD(anava, colnames(ad[i]), alpha=alpha.t)
                            letra1 <- letra$groups; colnames(letra1)=c("resp","groups")
                            if(transf !=1){letra1$respo=tapply(response,fatores[,i],mean, na.rm=TRUE)[rownames(letra1)]}}
                        print(letra1)
                        cat(green(bold("\n------------------------------------------")))
                        dadosm=data.frame(letra1,
                                          media=tapply(response, c(fatores[i]), mean, na.rm=TRUE)[rownames(letra1)],
                                          desvio=tapply(response, c(fatores[i]), sd, na.rm=TRUE)[rownames(letra1)])
                        dadosm$Tratamentos=factor(rownames(dadosm),levels = unique(unlist(fatores[i])))
                        dadosm$limite=dadosm$media+dadosm$desvio
                        dadosm=dadosm[as.character(unique(unlist(fatores[i]))),]
                        if(addmean==TRUE){dadosm$letra=paste(format(dadosm$media,digits = dec),dadosm$groups)}
                        if(addmean==FALSE){dadosm$letra=dadosm$groups}
                        media=dadosm$media
                        desvio=dadosm$desvio
                        Tratamentos=dadosm$Tratamentos
                        letra=dadosm$letra

                        grafico=ggplot(dadosm,aes(x=Tratamentos,
                                                  y=media))
                        if(fill=="trat"){grafico=grafico+
                            geom_col(aes(fill=Tratamentos),color=1)}
                        else{grafico=grafico+
                            geom_col(aes(fill=Tratamentos),fill=fill,color=1)}
                        if(errorbar==TRUE){grafico=grafico+
                            geom_text(aes(y=media+sup+if(sup<0){-desvio}else{desvio},
                                          label=letra),family=family,angle=angle.label, hjust=hjust)}
                        if(errorbar==FALSE){grafico=grafico+
                            geom_text(aes(y=media+sup,label=letra),family=family,angle=angle.label, hjust=hjust)}
                        if(errorbar==TRUE){grafico=grafico+
                            geom_errorbar(data=dadosm,
                                          aes(ymin=media-desvio,
                                              ymax=media+desvio,color=1),
                                          color="black",width=0.3)
                        grafico3=grafico+theme+
                            ylab(ylab)+
                            xlab(xlab)+
                            theme(text = element_text(size=textsize,color="black", family = family),
                                  axis.text = element_text(size=textsize,color="black", family = family),
                                  axis.title = element_text(size=textsize,color="black", family = family),
                                  legend.position = "none")
                        print(grafico3)}
                    }

                    #Para os fatores QUANTITATIVOS, regressao
                    if(quali[i]==FALSE && anavaF3[i,5]<=alpha.f){
                        cat(green(bold("\n------------------------------------------\n")))
                        cat('\nAnalyzing the simple effects of the factor ',fac.names[1],'\n')
                        cat(green(bold("\n------------------------------------------\n")))
                        cat(fac.names[i])
                        grafico3=polynomial(resp, fatores[,i])
                        cat(green("To edit graphical parameters, I suggest analyzing using the \"polynomial\" command"))
                    }

                    cat('\n')
                }
            }
        }

    #########################################################################################################################
    #Para interacao tripla significativa, desdobramento
    #########################################################################################################################
    if(anavaF3[8,5]<=alpha.f){
        cat(green(bold("\n------------------------------------------\n")))
        cat(green(bold("\nInteraction",paste(fac.names[1],'*',fac.names[2],'*',fac.names[3],sep='')," significant: unfolding the interaction\n")))
        cat(green(bold("\n------------------------------------------\n")))
        #===================================================================
        #Desdobramento de FATOR 1 dentro do niveis de FATOR 2 e do FATOR3
        #===================================================================

        cat(green(bold("\n------------------------------------------\n")))
        cat("Analyzing ", fac.names[1], ' within the combination of levels ', fac.names[2], 'and',fac.names[3])
        cat(green(bold("\n------------------------------------------\n")))

        m1=aov(resp~(Fator2*Fator3)/Fator1+bloco)
        anova(m1)
        pattern <- c(outer(levels(Fator2), levels(Fator3),
                           function(x,y) paste("Fator2",x,":Fator3",y,":",sep="")))
        des.tab <- sapply(pattern, simplify=FALSE,
                          grep, x=names(coef(m1)[m1$assign==5]))
        des1.tab <- summary(m1, split = list("Fator2:Fator3:Fator1" = des.tab))
        desd=des1.tab[[1]][-c(1,2,3,4,5),]
        desd=data.frame(desd[-length(rownames(desd)),])
        rownames(desd)=cbind(paste("Fator2:",rep(levels(Fator2),length(levels(Fator3))),
                                   "Fator3:",rep(levels(Fator3),e=length(levels(Fator2)))))
        colnames(desd)=c("Df",  "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
        print(desd)

        ii<-0
        for(i in 1:nv2) {
            for(j in 1:nv3) {
                ii<-ii+1
                # if(1-pf(QM/QME,glf,glE)[ii]<=alpha.f){
                if(quali[1]==TRUE){
                    cat('\n\n',fac.names[1],' inside of each level of ',lf2[i],' of ',fac.names[2],' and ',lf3[j],' of ',fac.names[3],"\n")
                    if(mcomp=='tukey'){tukey=TUKEY(resp[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],
                                                      fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],
                                                      anavaF3[9,1],
                                                      anavaF3[9,3],
                                                      alpha.t)
                    tukey=tukey$groups;colnames(tukey)=c("resp","letters")
                    if(transf !=1){tukey$respo=tapply(response[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],
                                                             fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],mean, na.rm=TRUE)[rownames(tukey)]}
                    # rownames(tukey)=levels(Fator1)[as.numeric(as.character(rownames(tukey)))]
                    print(tukey)}
                    if(mcomp=='duncan'){duncan=duncan(resp[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],
                                                      fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],
                                                      anavaF3[9,1],
                                                      anavaF3[9,3],
                                                      alpha.t)
                    duncan=duncan$groups;colnames(duncan)=c("resp","letters")
                    if(transf !=1){duncan$respo=tapply(response[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],
                                                              fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],mean, na.rm=TRUE)[rownames(duncan)]}
                    # rownames(tukey)=levels(Fator1)[as.numeric(as.character(rownames(tukey)))]
                    print(duncan)}
                    if(mcomp=='lsd'){lsd=LSD(resp[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],
                                                      fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],
                                                      anavaF3[9,1],
                                                      anavaF3[9,3],
                                                      alpha.t)
                    lsd=lsd$groups;colnames(lsd)=c("resp","letters")
                    if(transf !=1){lsd$respo=tapply(response[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],
                                                           fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],mean, na.rm=TRUE)[rownames(lsd)]}
                    # rownames(tukey)=levels(Fator1)[as.numeric(as.character(rownames(tukey)))]
                    print(lsd)}
                    if(mcomp=='sk'){
                        fat= fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]]
                        fat1=factor(fat,unique(fat))
                        levels(fat1)=1:length(levels(fat1))

                        sk=sk(resp[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],
                                    fat1,
                                    anavaF3$Df[9],
                                    anavaF3$`Sum Sq`[9],
                                    alpha.t)
                    colnames(sk)=c("resp","letters")
                    sk=sk[as.character(unique(fat1)),]
                    rownames(sk)=unique(fat)
                    if(transf !=1){sk$respo=tapply(response[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],
                                                   fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],mean, na.rm=TRUE)[rownames(sk)]}
                    # rownames(tukey)=levels(Fator1)[as.numeric(as.character(rownames(tukey)))]
                    print(sk)}
                    }
                # else{cat('\n\n',fac.names[1],' inside of each level of ',lf2[i],' of ',fac.names[2],' and ',lf3[j],' of ',fac.names[3])
                #     reg.poly(resp[fatores[,2]==lf2[i] & fatores[,3]==lf3[j]],
                #              fatores[,1][Fator2==lf2[i] & Fator3==lf3[j]],
                #              an[8,1],an[8,2], nv1-1, SQ[ii])}
                }
        }

        cat('\n\n')

        #===================================================================
        #Desdobramento de FATOR 2 dentro do niveis de FATOR 1 e FATOR 3
        #===================================================================

        cat("\n------------------------------------------\n")
        cat("Analyzing ", fac.names[2], ' within the combination of levels ', fac.names[1], 'and',fac.names[3])
        cat("\n------------------------------------------\n")
        m1=aov(resp~(Fator1*Fator3)/Fator2+bloco)
        anova(m1)
        pattern <- c(outer(levels(Fator1), levels(Fator3),
                           function(x,y) paste("Fator1",x,":Fator3",y,":",sep="")))
        des.tab <- sapply(pattern, simplify=FALSE,
                          grep, x=names(coef(m1)[m1$assign==5]))
        des1.tab <- summary(m1, split = list("Fator1:Fator3:Fator2" = des.tab))
        desd=des1.tab[[1]][-c(1,2,3,4,5),]
        desd=data.frame(desd[-length(rownames(desd)),])
        rownames(desd)=cbind(paste("Fator1:",rep(levels(Fator1),length(levels(Fator3))),
                                   "Fator3:",rep(levels(Fator3),e=length(levels(Fator1)))))
        colnames(desd)=c("Df",  "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
        print(desd)

        ii<-0
        for(k in 1:nv1) {
            for(j in 1:nv3) {
                ii<-ii+1
                #if(1-pf(QM/QME,glf,glE)[ii]<=alpha.f){
                if(quali[2]==TRUE){
                    cat('\n\n',fac.names[2],' inside of each level of ',lf1[k],' of ',fac.names[1],' and ',lf3[j],' of ',fac.names[3],'\n')
                    if(mcomp=='tukey'){tukey=TUKEY(resp[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],
                                                      fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]],
                                                      anavaF3[9,1],
                                                      anavaF3[9,3],
                                                      alpha.t)
                    tukey=tukey$groups;colnames(tukey)=c("resp","letters")
                    if(transf !=1){tukey$respo=tapply(response[fatores[,1]==lf1[i] & fatores[,3]==lf3[j]],
                                                             fatores[,2][Fator1==lf1[i]  & fatores[,3]==lf3[j]],mean, na.rm=TRUE)[rownames(tukey)]}
                    # rownames(tukey)=levels(Fator2)[as.numeric(as.character(rownames(tukey)))]
                    print(tukey)}
                    if(mcomp=='duncan'){duncan=duncan(resp[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],
                                                      fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]],
                                                      anavaF3[9,1],
                                                      anavaF3[9,3],
                                                      alpha.t)
                    duncan=duncan$groups;colnames(duncan)=c("resp","letters")
                    if(transf !=1){duncan$respo=tapply(response[fatores[,1]==lf1[i] & fatores[,3]==lf3[j]],
                                                              fatores[,2][Fator1==lf1[i]  & fatores[,3]==lf3[j]],mean, na.rm=TRUE)[rownames(duncan)]}

                    # rownames(tukey)=levels(Fator2)[as.numeric(as.character(rownames(tukey)))]
                    print(duncan)}
                    if(mcomp=='lsd'){lsd=LSD(resp[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],
                                                      fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]],
                                                      anavaF3[9,1],
                                                      anavaF3[9,3],
                                                      alpha.t)
                    lsd=lsd$groups;colnames(lsd)=c("resp","letters")
                    if(transf !=1){lsd$respo=tapply(response[fatores[,1]==lf1[i] & fatores[,3]==lf3[j]],
                                                           fatores[,2][Fator1==lf1[i]  & fatores[,3]==lf3[j]],mean, na.rm=TRUE)[rownames(lsd)]}
                    # rownames(tukey)=levels(Fator2)[as.numeric(as.character(rownames(tukey)))]
                    print(lsd)}
                    if(mcomp=='sk'){
                        fat=fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]]
                        fat1=factor(fat,unique(fat))
                        levels(fat1)=1:length(levels(fat1))
                        sk=sk(resp[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],
                                                  fat1,
                                                  anavaF3$`Sum Sq`[9],
                                                  anavaF3$`Sum Sq`[9],
                                                  alpha.t)
                    colnames(sk)=c("resp","letters")
                    sk=sk[as.character(unique(fat1)),]
                    rownames(sk)=unique(fat)
                    if(transf !=1){sk$respo=tapply(response[fatores[,1]==lf1[i] & fatores[,3]==lf3[j]],
                                                   fatores[,2][Fator1==lf1[i]  & fatores[,3]==lf3[j]],mean, na.rm=TRUE)[rownames(sk)]}
                    # rownames(tukey)=levels(Fator2)[as.numeric(as.character(rownames(tukey)))]
                    print(sk)}

                    }
                # else{cat('\n\n',fac.names[2],' inside of each level of ',lf1[k],' of ',fac.names[1],' and ',lf3[j],' of ',fac.names[3])
                #     reg.poly(resp[fatores[,1]==lf1[k] & fatores[,3]==lf3[j]],fatores[,2][Fator1==lf1[k] & fatores[,3]==lf3[j]], an[8,1], an[8,2], nv2-1, SQ[ii])}
                }
        }

        #===================================================================
        #Desdobramento de FATOR 3 dentro do niveis de FATOR 1 e FATOR 2
        #===================================================================

        cat(green(bold("\n------------------------------------------\n")))
        cat("Analyzing ", fac.names[3], ' within the combination of levels ', fac.names[1], 'and',fac.names[2])
        cat(green(bold("\n------------------------------------------\n")))

        m1=aov(resp~(Fator1*Fator2)/Fator3+bloco)
        anova(m1)
        pattern <- c(outer(levels(Fator1), levels(Fator2),
                           function(x,y) paste("Fator1",x,":Fator2",y,":",sep="")))
        des.tab <- sapply(pattern, simplify=FALSE,
                          grep, x=names(coef(m1)[m1$assign==5]))
        des1.tab <- summary(m1, split = list("Fator1:Fator2:Fator3" = des.tab))
        desd=des1.tab[[1]][-c(1,2,3,4,5),]
        desd=data.frame(desd[-length(rownames(desd)),])
        rownames(desd)=cbind(paste("Fator1:",rep(levels(Fator1),length(levels(Fator2))),
                                   "Fator2:",rep(levels(Fator2),e=length(levels(Fator1)))))
        colnames(desd)=c("Df",  "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
        print(desd)

        ii<-0
        for(k in 1:nv1) {
            for(i in 1:nv2) {
                ii<-ii+1
                # if(1-pf(QM/QME,glf,glE)[ii]<=alpha.f){
                if(quali[3]==TRUE){
                    cat('\n\n',fac.names[3],' inside of each level of ',lf1[k],' of ',fac.names[1],' and ',lf2[i],' of ',fac.names[2],'\n')
                    if(mcomp=='tukey'){tukey=TUKEY(resp[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                                      fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                                      anavaF3[9,1],
                                                      anavaF3[9,3],
                                                      alpha.t)
                    tukey=tukey$groups;colnames(tukey)=c("resp","letters")
                    if(transf !=1){tukey$respo=tapply(response[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                                             fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                                             mean, na.rm=TRUE)[rownames(tukey)]}
                    # rownames(tukey)=levels(Fator3)[as.numeric(as.character(rownames(tukey)))]
                    print(tukey)}
                    if(mcomp=='duncan'){duncan=duncan(resp[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                                      fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                                      anavaF3[9,1],
                                                      anavaF3[9,3],
                                                      alpha.t)
                    duncan=duncan$groups;colnames(duncan)=c("resp","letters")
                    if(transf !=1){duncan$respo=tapply(response[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                                              fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                                              mean, na.rm=TRUE)[rownames(duncan)]}
                    # rownames(tukey)=levels(Fator3)[as.numeric(as.character(rownames(tukey)))]
                    print(duncan)}
                    if(mcomp=='lsd'){lsd=LSD(resp[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                                      fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                                      anavaF3[9,1],
                                                      anavaF3[9,3],
                                                      alpha.t)
                    lsd=lsd$groups;colnames(lsd)=c("resp","letters")
                    if(transf !=1){lsd$respo=tapply(response[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                                           fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                                           mean, na.rm=TRUE)[rownames(lsd)]}
                    # rownames(tukey)=levels(Fator3)[as.numeric(as.character(rownames(tukey)))]
                    print(lsd)}
                    if(mcomp=='sk'){
                        fat=fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]]
                        fat1=factor(fat,unique(fat))
                        levels(fat1)=1:length(levels(fat1))
                        sk=sk(resp[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                                 fat1,
                                                 anavaF3$Df[9],
                                                 anavaF3$`Sum Sq`[9],
                                                 alpha.t)
                    colnames(sk)=c("resp","letters")
                    sk=sk[as.character(unique(fat1)),]
                    rownames(sk)=unique(fat)
                    if(transf !=1){sk$respo=tapply(response[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                                   fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                                                   mean, na.rm=TRUE)[rownames(sk)]}
                    print(sk)}

                    }
                # else{cat('\n\n',fac.names[3],' inside of each level of ',lf1[k],' of ',fac.names[1],' and ',lf2[i],' of ',fac.names[2])
                #     reg.poly(resp[fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                #              fatores[,3][fatores[,1]==lf1[k] & fatores[,2]==lf2[i]],
                #              an[8,1], an[8,2], nv3-1, SQ[ii])}
                }
        }

    }

    if(anavaF3[5,5]>alpha.f && anavaF3[6,5]>alpha.f && anavaF3[7,5]>alpha.f && anavaF3[8,5]>alpha.f){
        if(anavaF3[1,5]<=alpha.f | anavaF3[2,5]<=alpha.f | anavaF3[3,5]<=alpha.f){
        print(residplot)
        graficos}else{graficos=NA}}
    if(anavaF3[8,5]>alpha.f && anavaF3[5,5]<=alpha.f){
         graficos=list(residplot,colint1)
         if(anavaF3[6,5]>alpha.f && anavaF3[7,5]>alpha.f && anavaF3[3,5]<=alpha.f){
             graficos=list(residplot,colint1,grafico1)}
         graficos}
    if(anavaF3[8,5]>alpha.f && anavaF3[6,5]<=alpha.f){
         graficos=list(residplot,colint2)
         if(anavaF3[5,5]>alpha.f && anavaF3[7,5]>alpha.f && anavaF3[2,5]<=alpha.f){
             graficos=list(residplot,colint2,grafico2)}
         graficos}
    if(anavaF3[8,5]>alpha.f && anavaF3[7,5]<=alpha.f){
         graficos=list(residplot,colint3)
         if(anavaF3[5,5]>alpha.f && anavaF3[6,5]>alpha.f && anavaF3[1,5]<=alpha.f){
             graficos=list(residplot,colint3,grafico3)}}
    if(anavaF3[8,5]>alpha.f && anavaF3[5,5]<=alpha.f && anavaF3[6,5]<=alpha.f){
        graficos=list(residplot,colint1,colint2)
        graficos}
    if(anavaF3[8,5]>alpha.f && anavaF3[5,5]<=alpha.f && anavaF3[7,5]<=alpha.f){
        graficos=list(residplot,colint1,colint3)
        graficos}
    if(anavaF3[8,5]>alpha.f && anavaF3[6,5]<=alpha.f && anavaF3[7,5]<=alpha.f){
        graficos=list(residplot,colint2,colint3)
        graficos}
    if(anavaF3[8,5]<=alpha.f){graficos=list(residplot)}
    graficos=graficos
}

