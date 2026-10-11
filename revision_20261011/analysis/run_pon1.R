#!/usr/bin/env Rscript
# Bounded offline reanalysis; no credentials, network, manuscript edits or SNP selection.
options(stringsAsFactors=FALSE, digits=17)
base <- normalizePath(getwd())
out <- file.path(base, 'rerun_outputs', 'pon1')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
stopifnot(requireNamespace('TwoSampleMR', quietly=TRUE), requireNamespace('jsonlite', quietly=TRUE), requireNamespace('digest', quietly=TRUE))
sources <- c(olink=file.path(base, 'upstream_inputs', 'pon1', 'olink_PON1_chr7.gz'),
  finngen=file.path(base, 'upstream_inputs', 'pon1', 'finngen_R12_E4_AMYLNAS.gz'),
  api=file.path(base,'upstream_inputs','pon1','associations_direct.json'),
  provenance=file.path(base,'upstream_inputs','pon1','associations_direct.json.provenance.json'),
  metadata=file.path(base,'upstream_inputs','pon1','gwas_metadata.json'),
  script=file.path(base,'analysis','run_pon1.R'))
stopifnot(all(file.exists(sources)))
hash <- function(p) digest::digest(file=p,algo='sha256')
hashes <- vapply(sources,hash,character(1))
write.csv(data.frame(role=names(sources),path=unname(sources),sha256=unname(hashes)),file.path(out,'input_manifest.csv'),row.names=FALSE)
scan <- function(path,pattern,name,delimiter='') {
  con <- gzfile(path,'rt'); on.exit(close(con))
  header <- readLines(con,n=1,warn=FALSE); hits <- character(); total <- 0
  repeat { x <- readLines(con,n=100000,warn=FALSE); if(!length(x)) break
    total <- total+length(x); hits <- c(hits,x[grepl(pattern,x,fixed=TRUE)]) }
  writeLines(c(header,hits),file.path(out,paste0(name,'_source_rows.tsv')))
  cat(name,': scanned',total,'rows; matched',length(hits),'\n')
  stopifnot(length(hits)>0)
  # UKB-PPP is whitespace-delimited, whereas FinnGen is tab-delimited and may
  # contain spaces in text fields. Select the delimiter explicitly per source.
  read.table(text=paste(c(header,hits),collapse='\n'),header=TRUE,sep=delimiter,comment.char='',quote='',colClasses='character',check.names=FALSE)
}
p <- scan(sources[['olink']],'7:94937446:T:C:imp:v1','olink')
stopifnot(nrow(p)==1, p$ALLELE1=='C',p$ALLELE0=='T',p$GENPOS=='95308134',p$TEST=='ADD')
num <- function(x) as.numeric(x)
bx <- num(p$BETA); sx <- num(p$SE)
stopifnot(bx>0,sx>0,abs((bx/sx)^2-num(p$CHISQ))/num(p$CHISQ)<1e-5)
ex <- data.frame(SNP='rs662',beta.exposure=bx,se.exposure=sx,effect_allele.exposure='C',other_allele.exposure='T',eaf.exposure=num(p$A1FREQ),
  samplesize.exposure=num(p$N),exposure='PON1 UKB-PPP Olink discovery',id.exposure='Olink_PON1_rs662',mr_keep.exposure=TRUE)
# No underflowed P=0 is generated; source -log10(P) is retained separately.
write.csv(ex,file.path(out,'exposure_input.csv'),row.names=FALSE)
a <- jsonlite::fromJSON(sources[['api']]); a <- a[a$rsid=='rs662',]
ids <- c(T2D='ebi-a-GCST006867',HF='ebi-a-GCST009541',BMI='ieu-b-40')
stopifnot(nrow(a)==3,setequal(a$id,unname(ids)),!anyDuplicated(a$id),all(a$position==94937446),all(a$chr==7))
prov <- jsonlite::fromJSON(sources[['provenance']])
stopifnot(prov$sha256==hashes[['api']],prov$query$proxies==0)
o <- data.frame(SNP=a$rsid,beta.outcome=num(a$beta),se.outcome=num(a$se),effect_allele.outcome=a$ea,other_allele.outcome=a$nea,
  eaf.outcome=num(a$eaf),pval.outcome=num(a$p),id.outcome=a$id,outcome=names(ids)[match(a$id,ids)],mr_keep.outcome=TRUE)
f <- scan(sources[['finngen']],paste0('\t','rs662','\t'),'finngen',delimiter='\t')
f <- f[vapply(strsplit(f$rsids,'[,;]'),function(z) 'rs662' %in% z,logical(1)),]
stopifnot(nrow(f)==1,f[['#chrom']]=='7',f$pos=='95308134',setequal(c(f$ref,f$alt),c('C','T')))
o <- rbind(o,data.frame(SNP='rs662',beta.outcome=num(f$beta),se.outcome=num(f$sebeta),effect_allele.outcome=f$alt,other_allele.outcome=f$ref,
  eaf.outcome=num(f$af_alt),pval.outcome=num(f$pval),id.outcome='FinnGen_R12_E4_AMYLNAS',outcome='Amyloidosis',mr_keep.outcome=TRUE))
stopifnot(all(is.finite(o$beta.outcome)),all(o$se.outcome>0),all(o$effect_allele.outcome %in% c('C','T')))
write.csv(o,file.path(out,'outcome_inputs_before_harmonisation.csv'),row.names=FALSE)
h <- TwoSampleMR::harmonise_data(ex,o,action=2)
stopifnot(nrow(h)==4,all(h$mr_keep),!any(h$palindromic),all(h$effect_allele.outcome=='C'),all(h$other_allele.outcome=='T'))
write.csv(h,file.path(out,'harmonised_inputs.csv'),row.names=FALSE)
if(is.data.frame(attr(h,'log'))) write.csv(attr(h,'log'),file.path(out,'harmonisation_log.csv'),row.names=FALSE)
res <- TwoSampleMR::mr(h,method_list='mr_wald_ratio')
res <- res[match(c(ids,'FinnGen_R12_E4_AMYLNAS'),res$id.outcome),]
stopifnot(nrow(res)==4,all(res$nsnp==1))
res$CI95_lower <- res$b-qnorm(.975)*res$se
res$CI95_upper <- res$b+qnorm(.975)*res$se
res$P_Bonferroni_four_exploratory_tests <- p.adjust(res$pval,'bonferroni')
res$F_exposure <- (bx/sx)^2
res$exposure_N <- num(p$N)
res$exposure_scale <- 'Inverse-rank-normalised Olink protein level; not raw concentration'
res$outcome_scale <- ifelse(res$outcome=='BMI','GWAS-reported BMI scale','log odds')
res$interpretation <- 'Single-variant supplementary association; not independent cross-platform disease replication'
manual <- lapply(seq_len(nrow(res)),function(i) {
  z <- o[o$id.outcome==res$id.outcome[i],]
  sign <- if(z$effect_allele.outcome=='C' && z$other_allele.outcome=='T') 1 else if(z$effect_allele.outcome=='T' && z$other_allele.outcome=='C') -1 else stop('Unexpected alleles')
  b <- sign*z$beta.outcome/bx; s <- z$se.outcome/abs(bx); pv <- 2*pnorm(-abs(b/s))
  data.frame(id.outcome=z$id.outcome,original_effect_allele=z$effect_allele.outcome,flip=sign==-1,manual_b=b,manual_se=s,manual_p=pv,
    b_difference=b-res$b[i],se_difference=s-res$se[i],p_difference=pv-res$pval[i],source_gwas_P=z$pval.outcome,
    delta_SE_assuming_zero_covariance=sqrt(z$se.outcome^2/bx^2+z$beta.outcome^2*sx^2/bx^4))
})
checks <- do.call(rbind,manual)
stopifnot(all(abs(checks$b_difference)<1e-12),all(abs(checks$se_difference)<1e-12),all(abs(checks$p_difference)<1e-12))
# An additional orientation test: reversing both alleles and beta must preserve MR.
test_o <- o; test_o$effect_allele.outcome <- o$other_allele.outcome; test_o$other_allele.outcome <- o$effect_allele.outcome
test_o$beta.outcome <- -o$beta.outcome; test_o$eaf.outcome <- 1-o$eaf.outcome
test_h <- TwoSampleMR::harmonise_data(ex,test_o,action=2)
stopifnot(all(test_h$mr_keep),max(abs(test_h$beta.outcome-h$beta.outcome))<1e-12)
write.csv(res,file.path(out,'PON1_Wald_results.csv'),row.names=FALSE)
write.csv(checks,file.path(out,'independent_formula_checks.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
stopifnot(identical(hashes,vapply(sources,hash,character(1))))
jsonlite::write_json(list(run_time=as.character(Sys.time()),offline=TRUE,original_sources_unchanged=TRUE,all_four_retained=TRUE,
  manual_formula_checks_pass=TRUE,allele_reversal_test_pass=TRUE,OpenGWAS_snapshot=prov$retrieved_utc,
  framework_replication=FALSE,ITIH1_run=FALSE,manuscripts_modified=FALSE),file.path(out,'verification.json'),auto_unbox=TRUE,pretty=TRUE)
print(res[,c('outcome','b','se','CI95_lower','CI95_upper','pval','P_Bonferroni_four_exploratory_tests')])
