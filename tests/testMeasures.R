
## test linguistic measures

## data representation to the C function

##    location                     q1              q2                 q3
## 1:   Itziar          deitu,deitatu            diru freskatu,urreztatu
## 2:    Maule                  deitu             sos freskatu,urreztatu
## 3:  Senpere                  deitu          dihura             urritu
## 4:   Urketa deitatu,erran,atxikitu diru,sos,dihura   ihintzatu,urritu
## 5: Uztartze                               diru,sos                   
##             q4
## 1: herots,hots
## 2:      herots
## 3:            
## 4:            
## 5:   harrabots

set <- structure(list(gender = structure(c(1L, 2L, 2L, 1L, 2L, 1L, 1L, 
                          1L, 2L, 1L, 1L, 2L, 1L, 1L, 1L, 2L, 1L, 2L, 2L, 2L, 1L, 1L, 2L, 
                          1L, 1L, 2L), .Label = c("Female", "Male"), class = "factor"), 
                      location = c("Itziar", "Itziar", "Maule", "Urketa", "Urketa", 
                          "Urketa", "Senpere", "Itziar", "Maule", "Uztartze", "Uztartze", 
                          "Urketa", "Urketa", "Urketa", "Senpere", "Itziar", "Itziar", 
                          "Maule", "Maule", "Urketa", "Urketa", "Senpere", "Itziar", 
                          "Itziar", "Maule", "Uztartze"),
                      question = c("q1", "q1", "q1", "q1", "q1", "q1", "q1", "q2", "q2", "q2",
                          "q2", "q2", "q2", "q2", "q2", "q3", "q3", "q3", "q3", "q3", "q3", "q3", 
                          "q4", "q4", "q4", "q4"),
                      answer = c("deitu", "deitatu", "deitu", 
                          "deitatu", "erran", "atxikitu", "deitu", "diru", "sos", "diru", 
                          "sos", "diru", "sos", "dihura", "dihura", "freskatu", "urreztatu", 
                          "freskatu", "urreztatu", "ihintzatu", "urritu", "urritu", 
                          "herots", "hots", "herots", "harrabots")),
                 .Names = c("gender", "location", "question", "answer"),
                 row.names = c(NA, -26L), class = "data.frame")



ipiSimilarityMatrix <- matrix(
    c(0.683, 0.433, 0.5, 0.555, 
      0.666, 0.433, 0.5, 0.555,
      0.666, 0.3, 0.416, NA,
      0.35, 0.66, 0.166, NA,
      NA, 0.626, NA, 0.333),
    4, 5, TRUE, list(unique(set$question), unique(set$location)))


context('Measure correctness')

test_that('IRD correctness', {
              m <- diaMeasure(set, location ~ question, 'answer', 'ird', 'dice')
              m <- as.matrix(m)
              expect_that(round(m[1, 2], 2), equals(58.33)) # Itziar Maule (simple example)
              expect_that(round(m[4, 5], 2), equals(80))# Urketa Uztartze (NA observations)
              expect_that(round(m[4, 3], 2), equals(38.89)) # Sempere Urketa (more sample)
          })

test_that('IPD correctness', {
              set2 <- set[set$question %in% c("q1", "q1"), ]
              m2 <- diaMeasure(set2, location ~ question, 'answer', 'ipi', 'dice')
              m2 <- as.matrix(m2)
              expect_that(round(m[1, 2], 2), equals(67.41))
              
              m <- diaMeasure(set, location ~ question, 'answer', 'ipi', 'dice')
              m <- as.matrix(m)
          })
