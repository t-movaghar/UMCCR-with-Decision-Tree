# Time vector
t <- seq(0, 80, by = 0.1)

# --- Parameter Set 11: Label DT Untuned ---
c11 <- 0.28183867
alp1_11 <- 0.89473888
lam1_11 <- 0.22694818
alp2_11 <- 0.64982165
lam2_11 <- 0.01293865
alp3_11 <- 0.93945217
lam3_11 <- 0.02185218

H1_11 <- lam1_11 * t^alp1_11
H2_11 <- lam2_11 * t^alp2_11
H3_11 <- lam3_11 * t^alp3_11
S11 <- (c11 * exp(-H2_11 - H3_11)) + ((1 - c11) * exp(-H1_11 - H2_11 - H3_11))


# --- Parameter Set 1 dt tuned ---
c1 <- 0.28063302
alp1_1 <- 0.88079845
lam1_1 <- 0.21930927
alp2_1 <- 0.64093081
lam2_1 <- 0.01246792
alp3_1 <- 0.92455249
lam3_1 <- 0.02112281

H1_1 <- lam1_1 * t^alp1_1
H2_1 <- lam2_1 * t^alp2_1
H3_1 <- lam3_1 * t^alp3_1
S1 <- (c1 * exp(-H2_1 - H3_1)) + ((1 - c1) * exp(-H1_1 - H2_1 - H3_1))

# --- Parameter Set 2 rf tuned---
c2 <- 0.26497942
alp1_2 <- 0.87455577
lam1_2 <- 0.21531517
alp2_2 <- 0.63686662
lam2_2 <- 0.01222282
alp3_2 <- 0.91781909
lam3_2 <- 0.02075134

H1_2 <- lam1_2 * t^alp1_2
H2_2 <- lam2_2 * t^alp2_2
H3_2 <- lam3_2 * t^alp3_2
S2 <- (c2 * exp(-H2_2 - H3_2)) + ((1 - c2) * exp(-H1_2 - H2_2 - H3_2))

# --- Parameter Set 3 pals ---
c3 <- 0.204402264
alp1_3 <- 1.032766447
lam1_3 <- 0.273091477
alp2_3 <- 0.709014669
lam2_3 <- 0.013384210
alp3_3 <- 0.947013010
lam3_3 <- 0.017734819

H1_3 <- lam1_3 * t^alp1_3
H2_3 <- lam2_3 * t^alp2_3
H3_3 <- lam3_3 * t^alp3_3
S3 <- (c3 * exp(-H2_3 - H3_3)) + ((1 - c3) * exp(-H1_3 - H2_3 - H3_3))

# --- Parameter Set 4 svm tuned ---
c4 <- 0.25974165
alp1_4 <- 0.85941421
lam1_4 <- 0.20229590
alp2_4 <- 0.51070887
lam2_4 <- 0.02919325
alp3_4 <- 0.73210841
lam3_4 <- 0.04163230

H1_4 <- lam1_4 * t^alp1_4
H2_4 <- lam2_4 * t^alp2_4
H3_4 <- lam3_4 * t^alp3_4
S4 <- (c4 * exp(-H2_4 - H3_4)) + ((1 - c4) * exp(-H1_4 - H2_4 - H3_4))

# --- Plot ---
plot(t, S1, type = "l", col = "darkgreen", lwd = 2,
     xlab = "Time (in years)", ylab = "Overall survival probability",
     main = "Overall Survival Curves")
lines(t, S2, col = "blue", lwd = 2, lty = 2)
lines(t, S3, col = "purple", lwd = 2, lty = 3)
lines(t, S4, col = "orange", lwd = 2, lty = 4)
#lines(t, S5, col = "black", lwd = 2, lty = 5)

legend("topright",
       legend = c("DT Estimates", "RF Estimates", "Pal's Estimates", "SVM Estimates"),
       col = c("darkgreen", "blue", "purple", "orange"),
       lty = c(1, 2, 3, 4, 5),
       lwd = 2)


# --- Plot ---
plot(t, S3, type = "l", col = "purple", lwd = 2,
     xlab = "Time (in years)", ylab = "Overall survival probability",
     main = "Overall Survival Curves")
lines(t, S11, col = "darkgreen", lwd = 2, lty = 2)
#lines(t, S3, col = "purple", lwd = 2, lty = 3)

legend("topright",
       legend = c("Pal's Estimates", "DT Untuned Estimates"),
       col = c("purple", "darkgreen"),
       lty = c(1, 2, 3),
       lwd = 2)




# --- Parameter Set 6: DT with ADASYN ---
c6 <- 0.27661205
alp1_6 <- 0.87571583
lam1_6 <- 0.21818435
alp2_6 <- 0.67606960
lam2_6 <- 0.44213811
alp3_6 <- 0.91899207
lam3_6 <- 0.02102833

H1_6 <- lam1_6 * t^alp1_6
H2_6 <- lam2_6 * t^alp2_6
H3_6 <- lam3_6 * t^alp3_6
S6 <- (c6 * exp(-H2_6 - H3_6)) + ((1 - c6) * exp(-H1_6 - H2_6 - H3_6))

# --- Parameter Set 7: DT with Borderline ---
c7 <- 0.27595173
alp1_7 <- 0.87807474
lam1_7 <- 0.21822608
alp2_7 <- 0.39394391
lam2_7 <- 0.66233283
alp3_7 <- 0.92158810
lam3_7 <- 0.02102135

H1_7 <- lam1_7 * t^alp1_7
H2_7 <- lam2_7 * t^alp2_7
H3_7 <- lam3_7 * t^alp3_7
S7 <- (c7 * exp(-H2_7 - H3_7)) + ((1 - c7) * exp(-H1_7 - H2_7 - H3_7))

# --- Parameter Set 8: DT with SMOTE ---
c8 <- 0.26616351
alp1_8 <- 0.87375459
lam1_8 <- 0.21713696
alp2_8 <- 0.64712275
lam2_8 <- 0.46101063
alp3_8 <- 0.91689764
lam3_8 <- 0.02092491

H1_8 <- lam1_8 * t^alp1_8
H2_8 <- lam2_8 * t^alp2_8
H3_8 <- lam3_8 * t^alp3_8
S8 <- (c8 * exp(-H2_8 - H3_8)) + ((1 - c8) * exp(-H1_8 - H2_8 - H3_8))

# --- Parameter Set 9: DT with Random Undersampling ---
c9 <- 0.29629630
alp1_9 <- 0.61445574
lam1_9 <- 0.02542384
alp2_9 <- 0.50215815
lam2_9 <- 0.03150710
alp3_9 <- 0.81181273
lam3_9 <- 0.01737444

H1_9 <- lam1_9 * t^alp1_9
H2_9 <- lam2_9 * t^alp2_9
H3_9 <- lam3_9 * t^alp3_9
S9 <- (c9 * exp(-H2_9 - H3_9)) + ((1 - c9) * exp(-H1_9 - H2_9 - H3_9))

# --- Parameter Set 10: DT with Random Oversampling ---
c10 <- 0.2576124
alp1_10 <- 0.8773073
lam1_10 <- 0.2162088
alp2_10 <- 0.6625802
lam2_10 <- 0.3008548
alp3_10 <- 0.9055170
lam3_10 <- 0.2067561

H1_10 <- lam1_10 * t^alp1_10
H2_10 <- lam2_10 * t^alp2_10
H3_10 <- lam3_10 * t^alp3_10
S10 <- (c10 * exp(-H2_10 - H3_10)) + ((1 - c10) * exp(-H1_10 - H2_10 - H3_10))


# --- Plot ---
plot(t, S1, type = "l", col = "darkgreen", lwd = 2,
     xlab = "Time (in years)", ylab = "Overall survival probability",
     main = "Survival Curves for Decision Trees with Resampling Methods")
  ylim = c(0, 1)

lines(t, S6, col = "blue", lwd = 2, lty = 2)       # ADASYN
lines(t, S7, col = "purple", lwd = 2, lty = 3)     # Borderline-SMOTE
lines(t, S8, col = "orange", lwd = 2, lty = 4)     # SMOTE
lines(t, S9, col = "red", lwd = 2, lty = 5)        # Undersampling
lines(t, S10, col = "green", lwd = 2, lty = 6)     # Oversampling
lines(t, S3, col = "yellow", lwd = 2, lty = 7)     # Oversampling

legend("topright",
       legend = c("DT None", "DT ADASYN", "DT Borderline", "DT SMOTE", "DT Undersampling", "DT Oversampling", "Pal's Estimates"),
       col = c("darkgreen", "blue", "purple", "orange", "red", "green", "yellow"),
       lty = 1:7,
       lwd = 2)
