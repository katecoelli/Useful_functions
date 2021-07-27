#####################################################################
#################### Soil related functions #########################
#####################################################################

# This code is a combination of functions written by various people #
# If I don't know the original author I have attributed it to the person
# who gave me the code


### GLOBAL BD PTF ###

# This PTF is the PTF by Tranter et al, 2007
Van_Bemelen_factor<- 1.72

#function to calculate bulk density - global PTF
bd_glob <- function(OC, Van_Bemelen_factor, sand, mid_depth){
  om=OC*Van_Bemelen_factor
  100 /(om/0.223  + (100 - om) /
          (1.35128606477631 + 0.00451974677070142 * sand + 
             (sand - 44.652494600432) * ((sand - 44.652494600432) *
                                           -0.0000613723924995459) + 0.0596420803366252 * log(mid_depth)))
}  



### FIELD CAPACITY ###

#function was given to me from James Maloney - originally derived in 
#Jose Campusano's masters thesis

FC<- function(sand,clay){
  0.4795 - 3.873 * 10^-5 * sand ^2 - 6.701 * 10^-7 * clay ^2 * sand
}



### PERMANENT WILTING POINT ###

#function was given to me from James Maloney - originally derived in 
#Jose Campusano's masters thesis

PWP<- function(sand, clay, mid_depth, BD){
  -0.1554 - 0.7221 * tanh(0.5 * (-0.9705 - 0.8529 * BD - 0.00827 *
                                   clay + 0.01994 * sand))  + 0.1325 * tanh(0.5 * (3.71 - 3.19 * BD + 0.01205 * clay + 0.01617 * sand)) + 
    0.1720 * tanh(0.5 * (-3.94 - 0.5067 * BD + 0.02158 * clay + 0.04978 * sand))
  
}

