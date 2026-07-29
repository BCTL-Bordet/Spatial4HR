## forest plot

formatNdig = function(x, n=2)
{ m = floor(log10(x))+1;
if (m>=n) { return(format(x, digits=2)); }  ########### CHANGED from digits = 0 to = 2 due to error: Error in prettyNum(.Internal(format(x, trim, digits, nsmall, width, 3L, invalid value 0 for 'digits' argument
#m2 = m; if (m2<0) { m2 = 0; }
format(x, nsmall=1, digits=n) ### CHANGED from nsmall=n-m to 1 for same error above
}

formatNice = function(x, parse=TRUE, lim=1e-20, nDigits=2)
{ if (!is.finite(x))
{ if (is.na(x)) { return("NA"); }
  if (x==Inf) { return("Inf"); }
  if (x==-Inf) { return("-Inf"); }
  return("NA");
}
  if (x < lim) { return("0"); }
  if (abs(x)>=1e-3 && abs(x)<=1e3) { return(formatNdig(x, n=nDigits)); }
  ex = ceiling(-log10(x));
  ma = formatNdig(x*10^ex, n=nDigits);
  #ret = paste(ma,"%*%10^{",-ex,"}");
  if (ma=="1") { ret=paste0("10^{",-ex,"}")}
  else { ret = paste0("paste('", ma, "x', 10^{",-ex,"})"); }
  if (parse) { return(parse(text=ret)); } 
  return(ret);
}


formatNiceP = function(x, descr=NULL, parse=TRUE, lim=1e-20, pName="p", sepDescr=" - ", boldDescr=FALSE)
{ p = formatNice(x, parse=FALSE, lim=lim);
if (!grepl("\\^",p)) { p = paste0("'", p, "'")}
if (is.null(descr)) { ret = paste0("paste(italic(", pName, "),' = ', ", p, ")"); }
else
{ if (boldDescr) { descr = paste0("bold('", descr, "')"); } else { descr = paste0("'", descr, "'")}
  ret = paste0("paste(", descr, ",'", sepDescr, "', italic(", pName, "),' = ', ", p, ")");
}
if (parse) { return(parse(text=ret)); } 
return(ret);
}


fpDrawCI = function (lower_limit, estimate, upper_limit, size, y.offset = 0.5, 
                     clr.line, clr.marker, lwd = 1.5, lty = 1, vertices, vertices.height = 0.1, shapes_gp = fpShapesGp(), 
                     shape_coordinates = structure(c(1, 1), max.coords = c(1, 1)),
                     ...) 
{
  forestplot:::prFpDrawLine(lower_limit = lower_limit, upper_limit = upper_limit, 
                            clr.line = clr.line, lwd = lwd, lty = lty, y.offset = y.offset, 
                            vertices = vertices, vertices.height = vertices.height)
  box <- convertX(unit(estimate, "native"), "npc", valueOnly = TRUE)
  skipbox <- box < 0 || box > 1
  if (!skipbox) {
    if (!is.unit(size)) {
      size <- unit(size, "snpc")
    }
    # grid.rect(x = unit(estimate, "native"), y = y.offset,                                                        # this will plot squares
    #     width = size, height = size, gp = gpar(fill = clr.marker, 
    #         col = clr.line));
    grid.circle(x = unit(estimate, "native"), y = y.offset, r = size,                                              # this will plot circles instead of squares
                gp = prGetShapeGp(shapes_gp, shape_coordinates, "box", default = gpar(fill = clr.marker,  
                                                                                      col = clr.marker)))
  }
}


# updated with names_var and arrows_surv_pCR options
allForest = function(x, y, w=NULL, control=~1, sig=.05, sigOnFDR=FALSE, fdr=TRUE, new_page=FALSE, clip=c(.24,4.1),        # sigOnFDR = TRUE will highlight only things with FDR significant
                     ticks=c(.25,.5,1,2,4), dispN=FALSE, useWilcox=FALSE, genesItal=FALSE, parseNames=FALSE, rmNAp=TRUE, normReal=TRUE, 
                     direction_colors = TRUE,                                         # direction_colors=TRUE will (should!) give red=worse (HR>1 or OR<1) and blue = better (HR<1 or OR>1)
                     boxsize = 0.1,                                                      # boxsize = 0.3 (for example) will fix the size of the circles/squares
                     graphwidth = unit(4, "cm"),                                        # graphwidth = unit(4, "cm") (for example) will fix the size of the plot (but not the "table" on the left with labels/HR/CI/P value/FDR!)
                     capHR=TRUE, 
                     names_var= NA, # we can give new names to the variables in plot 
                     arrows_surv_pCR = TRUE, # if working with survival analyses or pCR leave TRUE, it adds the arrows in the bottom part of the plot. If FALSE, it just adds HR or OR
                     ...)    
{ #f = function(x) format(x, digit=2);
  if (!require(forestplot)) { stop("Need forestplot"); }
  f = formatNice;
  plBlue = function(lower_limit, estimate, upper_limit, size, y.offset=.5, ...)
  { fpDrawCI(lower_limit, estimate, upper_limit, size, y.offset=.5, 
             # clr.line="darkblue", 
             clr.line="#56B4E9",       # "#00AFBB" , "royalblue1                ################# In this three functions (plBlue, plYellow, plRed) you can change the colors.
             clr.marker="#56B4E9"
             ,vertices = TRUE                  # added
  );
  }
  plYellow = function(lower_limit, estimate, upper_limit, size, y.offset=.5, ...) # not yellow anymore, but same name of the function
  { fpDrawCI(lower_limit, estimate, upper_limit, size, y.offset=.5, 
             clr.marker = alpha("slategray4", 1), # "orange" , "darkgrey"
             clr.line="slategray4"
             # clr.line="yellow4" 
             ,vertices = TRUE
  );
  }
  plRed = function(lower_limit, estimate, upper_limit, size, y.offset=.5, ...)
  { fpDrawCI(lower_limit, estimate, upper_limit, size, y.offset=.5, 
             # clr.line="darkblue", 
             clr.line="#D55E00",        # "#FC4E07", "red4"
             clr.marker="#D55E00"
             ,vertices = TRUE
  );
  }
  if (class(y) == "Surv") { xl = "HR"; } else { xl = "OR"; }
  
  if (arrows_surv_pCR) { 
    if (class(y) == "Surv") { xarrow = expression("" [italic("Better outcome ") %<-% "  "]* bold("HR") ["  " %->% italic(" Worse outcome")]); } 
    else { xarrow = expression("" [italic("Less pCR ") %<-% "  "]* bold("OR") ["  " %->% italic(" More pCR")]); }
  }
  else{ 
    if (class(y) == "Surv") { xarrow = expression(bold("HR")); } else { xarrow = expression(bold("OR")); } 
  }
  
  
  fbO = update(y~1, control);
  
  nms = colnames(x);
  if (is.data.frame(x) & normReal)
  { x = lapply(x, function(i) if (is.numeric(i) & !is.integer(i))
  { (i-mean(i, na.rm=TRUE))/sd(i, na.rm=TRUE); } else { i; } );
  x = do.call(data.frame, x);
  }
  else { x = data.frame(x); }
  names(nms) = names(x);
  
  if (is.null(w)) { w = rep(TRUE, nrow(x)); }
  w = w & !is.na(y);
  if (length(all.vars(control))>0)
  { w = w & !rowAnys(is.na(x[,all.vars(control),drop=FALSE]));
  }
  x = x[w,];
  y = y[w];
  
  vOk = setdiff(colnames(x), all.vars(control));
  
  tt = list();
  tt[[1]] = as.list(c("", nms[vOk]));
  if (genesItal)
  { if (!require(org.Hs.eg.db)) { stop("Need library org.Hs.eg.db"); }
    w = tt[[1]] %in% keys(org.Hs.eg.db, "ALIAS");
    for (i in which(w)) { tt[[1]][[i]] = parse(text=paste0("italic('",tt[[1]][[i]],"')")) }
  }
  if (parseNames)
  { tt[[1]] = lapply(tt[[1]], function(i) parse(text=i))
  }
  vls = list("", parse(text=paste0("bold(",xl,")")), expression(bold("95% CI")),
             expression(bolditalic(p)), expression(bold(FDR)));
  if (!fdr) { vls=vls[-length(vls)]; }
  for (i in 2:length(vls))
  { tt[[i]] = list(vls[[i]]);
  }
  
  tn = list();
  col = list();#rep("black", length(vOk));
  pp = rep(NA, length(vOk));
  for (i in seq_along(vOk))
  { fb = fbO;
  if (!identical(control, ~1))
  { if (any(is.na(x[,vOk[i]]))) # Look if some parameters got collapsed
  { av = all.vars(fb)[-1];
  Noks = sapply(av, function(id) length(unique(x[[id]][!is.na(x[,vOk[i]])])))<2;
  for (id in names(which(Noks))) { fb = update(fb, formula(paste("~ . -", id))); }
  }
    if (class(y) == "Surv")
    { a0 = coxph(fb, data=x, subset=!is.na(x[,vOk[i]])); }
    else { a0 = glm(fb, data=x, family=binomial, subset=!is.na(x[,vOk[i]])); }
  }
  
  fn = update(fb, formula(paste("~", vOk[i],"+.")));
  if (class(y) == "Surv")
  { a = coxph(fn, data=x);
  if (identical(control, ~1))
  { if (useWilcox && (is.factor(x[,vOk[i]]) || is.logical(x[,vOk[i]])))
  { fm = survdiff(fn, data=x)
  p = pchisq(fm$chisq, length(fm$n) - 1, lower.tail = FALSE)
  }
    else { p = summary(a)$logtest["pvalue"]; }
  }
  else { an = anova(a0, a); p = an[["Pr(>|Chi|)"]][2]; }  ### changed with R 4.2.1 from P(>|Chi|) to Pr(>|Chi|)
  }
  else
  { a = glm(fn, data=x, family=binomial);
  if (identical(control, ~1))
  { if (useWilcox)
  { if (is.numeric(x[,vOk[i]])) { p = wilcox.test(x[unclass(factor(y))==1,vOk[i]], x[unclass(factor(y))==2,vOk[i]])$p.value; }
    else { p = fisher.test(x[,vOk[i]], y)$p.value; }
  }
    else
    { p = coef(summary(a));
    p = p[nrow(p),"Pr(>|z|)"];
    }
  }
  else
  { an = anova(a0, a);
  p = pchisq(an$Deviance[2], an$Df[2], lower.tail=FALSE);
  }
  }
  
  co = coef(summary(a));
  if (class(y) == "Surv") { me = co[1,"coef"]; se = co[1,"se(coef)"]; }
  else { me = co[2,1]; se = co[2,2]; }
  if (!rmNAp || !is.na(me))
  { if (useWilcox && !class(y) == "Surv" && !is.numeric(x[,vOk[i]]) && length(setdiff(unique(x[,vOk[i]]), NA))==2)
  { fm = fisher.test(x[,vOk[i]], y);
  tt[[2]][[i+1]] = f(fm$estimate);
  tt[[3]][[i+1]] = parse(text=paste("paste(",f(fm$conf.int[1]), ",' to ',", f(fm$conf.int[2]), ")"));
  tn[[i]] = log(c(fm$estimate, fm$conf.int));
  }
    else
    { if (capHR & abs(me)>5) {  tt[[2]][[i+1]] = ifelse(me<0, 0, "+Inf"); me = ifelse(me<0, log(clip[1]), log(clip[2])); }
      else { tt[[2]][[i+1]] = f(exp(me)); }
      #tt[[3]][[i+1]] = paste(f(exp(me-1.96*se)), "to", f(exp(me+1.96*se)));
      tt[[3]][[i+1]] = parse(text=paste("paste(",f(exp(me-1.96*se)), ",' to ',", f(exp(me+1.96*se)), ")"));
      tn[[i]] = c(me, me-1.96*se, me+1.96*se);
    }
    tt[[4]][[i+1]] = f(p);
    pp[i] = p;
  } else { tn[[i]] = c(NA, NA, NA); pp[i] =NA; tt[[2]][[i+1]]=tt[[3]][[i+1]]=tt[[4]][[i+1]]=""; }
  }
  tn = do.call(rbind, tn);
  colnames(tn) = c("mean", "lower", "upper");
  tn = exp(tn);
  
  tn = rbind(NA, tn);
  
  if (fdr)
  { fdr = p.adjust(pp, method='fdr');
  tt[[5]][2:(length(fdr)+1)] = lapply(fdr,function(i) if (rmNAp & is.na(i)) { "" } else f(i));
  }
  
  
  if (!is.null(sig)) {
    if (direction_colors) {
      if (class(y) == "Surv") {
        if (sigOnFDR) { psig = fdr; } else { psig = pp; }
        dir_col = as.numeric(tt[[2]][2:length(tt[[2]])])
        col = mapply(psig, dir_col, FUN = function(i, z) if (!is.na(i) & i < sig & z < 1) { plBlue } else if (!is.na(i) & i < sig & z > 1) { plRed } else { plYellow });
        
      }
      else { 
        if (sigOnFDR) { psig = fdr; } else { psig = pp; }
        dir_col = as.numeric(tt[[2]][2:length(tt[[2]])])
        col = mapply(psig, dir_col, FUN = function(i, z) if (!is.na(i) & i < sig & z > 1) { plBlue } else if (!is.na(i) & i < sig & z < 1) { plRed } else { plYellow });
      }
    }
    
    else  {
      if (sigOnFDR) { psig = fdr; } else { psig = pp; }
      col = lapply(psig, function(i) if (!is.na(i) & i < sig) { plBlue } else { plYellow });
    }
  }
  
  if (dispN)
  { N = colSums(!is.na(x[,vOk]));
  tt = c(tt[1], list(c(list(expression(bold(N))), as.list(N))), tt[2:length(tt)]);
  }
  
  is = c(FALSE, rep(FALSE, ncol(x)));
  col = c(col[[1]], col);
  
  wn = which(tn[,3]==tn[,2]);
  tn[wn,] = NA;
  
  attr(ticks, "labels")=as.character(ticks);
  if(!identical(names_var, NA)){
    tt[[1]] = as.list(c("", names_var));}
  else {tt[[1]] = as.list(tt[[1]])
  }
  print(grob<-forestplot(tt, tn, new_page = new_page, xlog=TRUE, is.summary=is, xticks=ticks,
                         fn.ci_norm=col, 
                         
                         # xlab=xl,   # original, not used anymore
                         xlab=xarrow, # if arrows_surv_pCR = TRUE working with survival analyses or pCR (it adds the arrows in the bottom part of the plot)
                         
                         clip=clip,
                         col = fpColors(zero = "black"),
                         # txt_gp=fpTxtGp(xlab=gpar(cex=.7), ticks=gpar(cex=.7)), ...));
                         txt_gp=fpTxtGp(xlab=gpar(cex=1), ticks=gpar(cex=1), label = gpar(cex = 1)),    # change this if you want different size of labels/text in general
                         boxsize = boxsize,
                         graphwidth = graphwidth
                         , ...));
  return(invisible(list(tt=tt, tn=tn, grob=grob)));
}

