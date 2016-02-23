library(GEOquery)
gse <- 'GSE50192'
geo <- getGEO(gse, GSEMatrix=T)
raw.geo <- getGEOSuppFiles(gse)
geo2 <- readGEORawFile(file.path(gse, 'GSE50192_Matrix_Signal.txt'), sep="\t", Uname='Signal_A', Mname='Signal_B', mergeManifest=T)
m.geo <- getM(geo2)
mval <- m.geo
mval <- mval[!apply(is.na(mval), 1, any),]
design <- model.matrix(int.formula, data=targets2)

ctrls <- getProbeInfo(raw.meth, type='Control')
cg.ctrls <- ctrls$Address %in% array.annot$AddressA

inc <- getINCs(geo2)
m.inc <- rbind(mval, inc)
ctl <- rownames(m.inc) %in% rownames(inc)
fit1 <- RUVfit(data=m.inc, design=design, coef=2, ctl=ctl, method=method)
fit2 <- RUVadj(fit1)
top1 <- topRUV(fit2, number=Inf)
ctl <- rownames(mval) %in% rownames(top1[top1$p.ebayes.BH > 0.5,])
fit1 <- RUVfit(data=mval, design=design, coef=2, ctl=ctl, method=method)
fit2 <- RUVadj(fit1)
top2 <- topRUV(fit2, number=Inf)
fitvar <- varFit(mval, design=design, coef=c(1,2))
topDV <- topVar(fitvar, coef=2, num=nrow(top2))
mval.clean <- removeBatchEffect(mval, covariates=t(fit1$W), design=design)
beta.clean <- 2^mval.clean/(1+2^mval.clean)



