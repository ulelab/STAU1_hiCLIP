hotfactor= function(df,n=10,o="other") {
  # copied from : http://stackoverflow.com/questions/6673574/how-can-i-replace-a-factor-levels-with-the-top-n-levels-by-some-metric-plus
   fac <- df$Var1
   by <- df$Freq
   levels(fac)[rank(-xtabs(by~fac))[levels(fac)]>n] <- o
   fac <- reorder(fac,-by,sum)
   fac.level <- levels(fac)
   fac <- factor(fac, levels = c(as.character(fac.level[!(fac.level %in% o)]), o))
   fac
}