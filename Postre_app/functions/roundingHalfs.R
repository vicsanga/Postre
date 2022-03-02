##To improve rounding limitations, because .5 are rounded to closer even number... so 1.5 rounded to 2 and 0.5 rounded to 0
##That is why janitor package implements round_half_up function
##+info
# https://stackoverflow.com/questions/12688717/round-up-from-5/12688836#12688836
##https://cran.r-project.org/web/packages/janitor/janitor.pdf


# In base R round(), halves are rounded to even, e.g., 12.5 and 11.5 are both rounded to 12.  Thisfunction rounds 12.5 to 13 (assumingdigits = 0).  
# Negative halves are rounded away from zero,e.g., -0.5 is rounded to -1.This may skew subsequent statistical analysis of the data, but may be desirable in certain contexts.
# This function is implemented exactly fromhttps://stackoverflow.com/a/12688836; see thatquestion and comments for discussion of this issue

##The solution in stack overflow
#This is the code copy pasted from the janitor package round_half_up function
#Same signature that R built in round() function
round2 = function (x, digits = 0) {
  posneg <- sign(x)
  z <- abs(x) * 10^digits
  z <- z + 0.5 + sqrt(.Machine$double.eps)
  z <- trunc(z)
  z <- z/10^digits
  z * posneg
}



