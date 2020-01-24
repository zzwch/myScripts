install.packages("Vennerable", repos="http://R-Forge.R-project.org")

library(Vennerable)

CD44 = 1:1000
PK44= 1:111
NE = c(27:111, 200:557)


venn <- Venn(Sets = list(CD44 = CD44,
                         NE = NE,PK44 = PK44))
plot(venn[,c("NE", "PK44")],
     doWeight = T)
plot(compute.Venn(venn[,c("NE", "PK44")], 
                  doWeights = TRUE, 
                  doEuler = T,
                  type = "circles"),
     show = list(
       #FaceText = "signature",
       SetLabels = T,
       Faces = F,
       DarkMatter = T
     )
)

