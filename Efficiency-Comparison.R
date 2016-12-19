#YU_WANG_Final_Project

#1
dose<- read.csv("doseresponsedata16.csv")
rfinal<- dose[dose$RFinal != "*",]
rf<- as.numeric(levels(rfinal$RFinal))[rfinal$RFinal]
r1<- as.numeric(levels(rfinal$R1))[rfinal$R1]
r2<- as.numeric(levels(rfinal$R2))[rfinal$R2]
r3<- as.numeric(levels(rfinal$R3))[rfinal$R3]
rfinalA<- as.data.frame(matrix(c(rfinal$Patient,rfinal$Dose,rf),ncol = 3))
colnames(rfinalA) <- c("Patient","Dose","rf")
dose0<-(rfinalA[rfinal$Dose == "0",])
dose2<-(rfinalA[rfinal$Dose == "2",])
dose10<-(rfinalA[rfinal$Dose == "10",])
dose50<-(rfinalA[rfinal$Dose == "50",])

#outlier check
boxplot(rfinalA$rf ~ rfinalA$Dose,main="By Dose Level")
boxplot(rfinalA$rf,main="Overall RFinal")
install.packages("outliers")
library(outliers)
chisq.out.test(dose10$rf, variance=var(dose10$rf), opposite = FALSE)
chisq.out.test(rfinalA$rf, variance=var(rfinalA$rf), opposite = FALSE)

rfinalA<- rfinalA[rfinalA$rf != 8.5,]
dose10<- dose10[dose10$rf != 8.5,]

#normality check
par(mfrow=c(1,2))
hist(rf,main="histogram of RFinal",xlab="RFnial")
qqnorm(rf,main="QQ Plot of RFinal")
qqline(rf)
par(mfrow=c(2,2))
hist(dose0$rf,main="histogram of dose0 RFinal",xlab="RFnial")
hist(dose2$rf,main="histogram of dose2 RFinal",xlab="RFnial")
hist(dose10$rf,main="histogram of dose10 RFinal",xlab="RFnial")
hist(dose50$rf,main="histogram of dose50 RFinal",xlab="RFnial")

qqnorm(dose0$rf,main="QQ Plot of dose0 RFinal")
qqline(dose0$rf)
qqnorm(dose2$rf,main="QQ Plot of dose2 RFinal")
qqline(dose2$rf)
qqnorm(dose10$rf,main="QQ Plot of dose10 RFinal")
qqline(dose10$rf)
qqnorm(dose50$rf,main="QQ Plot of dose50 RFinal")
qqline(dose50$rf)
shapiro.test(rfinalA$rf)
shapiro.test(dose0$rf)
shapiro.test(dose2$rf)
shapiro.test(dose10$rf)
shapiro.test(dose50$rf)

#variance
install.packages("lawstat")
library(lawstat)
bartlett.test(rfinalA$rf ~ rfinalA$Dose)
levene.test(rfinalA$rf,rfinalA$Dose)

#(a)
#contrast
#ANOVA test
rfinalA$Dose = factor(rfinalA$Dose, levels=unique(rfinalA$Dose))
levels(rfinalA$Dose)
vs1= c(-4,2,6,-4)
vs2= c(-1,-1,3,-1)
ctr <- cbind(vs1,vs2)
contrasts(rfinalA$Dose) = ctr
CList = list("vs1" = 1,"vs2" = 2)
model = aov(rf ~ Dose, data = rfinalA)             
summary(model, split=list(Dose=CList))
aov.out = aov(rf ~ Dose, data = rfinalA)
summary(aov.out)
TukeyHSD(aov.out)
vsd<- TukeyHSD(aov.out)$Dose[c(2,5,6),]
vsd[1,]<- -vsd[1,]
vs<- rbind(-4*vsd[1,], -2*vsd[2,],-2*vsd[3,])
colSums(vs)

#t test 
rfinalA$Dose = factor(rfinalA$Dose, levels=unique(rfinalA$Dose))
levels(rfinalA$Dose)
m1 = lm(rf ~ Dose,data = rfinalA)
Input =
  "c              Dose2  Dose10 Dose50 Dose0
    dvs1            -4     2   6   -4   
    dvs2             -1    -1   3   -1    
    "
Matriz = as.matrix(read.table(textConnection(Input), header=TRUE,row.names=1))
library(multcomp)
t <- glht(m1, linfct = Matriz)
t$linfct
summary(t)
vsd1<- TukeyHSD(aov.out)$Dose[c(2,4,6),]
vsd1[3,]<- -vsd1[3,]
colSums(vsd1)

#(b)
rfinalA$Dose = factor(rfinalA$Dose)
rfinalA$Dose <- relevel(rfinalA$Dose, "0")
trend <- lm(rf ~ Dose, data = rfinalA)
library(multcomp)
library(sandwich)
sandwich(trend)
covmatrix <- vcovHC(trend, type = "HC") # Heteroskedasticity-Consistent Covariance Matrix Estimation
Williams.trend <- glht(trend, linfct=mcp(Dose = "Williams"),vcov=sandwich, alternative="less")
par(mfrow=c(1,1))
plot(Williams.trend, main="Williams trend approach for decreasing trend")
summary(Williams.trend)

#(c)
n <- tapply(rfinalA$Dose, rfinalA$Dose, length)
k <- length(n)
Contrast.M <- c()
for (i in 1:(k - 1)) {
  help <- c(-1, n[2:(i + 1)] / sum(n[2:(i + 1)]), rep(0 , k - i - 1))
  Contrast.M <- rbind(Contrast.M, help)
}
rownames(Contrast.M) <- paste("C", 1:nrow(Contrast.M))
Contrast.M
Tmax <- glht(trend, linfct=mcp(Dose = Contrast.M),alternative="less")
summary(Tmax)

#(d)
r.matrix<-cbind(r1,r2,r3,rf)
rf.mean<-rowSums (r.matrix)/4
rf.d<- cbind(rfinal,rf.mean)
rf.d<-rf.d[,c(1,2,3,4,9)]


rf.d$Dose = factor(rf.d$Dose)
levels(rf.d$Dose)
trend.d <- lm(rf.mean ~ Dose+Bvalue1+Bvalue2, data = rf.d)
library(multcomp)
library(sandwich)
sandwich(trend.d)
covmatrix.d <- vcovHC(trend.d, type = "HC") # Heteroskedasticity-Consistent Covariance Matrix Estimation
Williams.trend.d <- glht(trend.d, linfct=mcp(Dose = "Williams"),vcov=sandwich, alternative="less")
plot(Williams.trend.d, main="Williams trend approach for decreasing trend with covariate included")
summary(Williams.trend.d)

n <- tapply(rf.d$Dose, rf.d$Dose, length)
k <- length(n)
CM <- c()
for (i in 1:(k - 1)) {
  help <- c(-1, n[2:(i + 1)] / sum(n[2:(i + 1)]), rep(0 , k - i - 1))
  CM <- rbind(CM, help)
}
rownames(CM) <- paste("C", 1:nrow(CM))
CM
Tmax.d <- glht(trend.d, linfct=mcp(Dose = CM),alternative="less")
summary(Tmax.d)

#2
app1<- read.csv("MAT750.app1.csv")

#(a)
#data summary
par(mfrow=c(2,5))
plot(app1$gender,app1$app_1,main="gender vs app1_usage")
plot(app1$gender,app1$app_2,main="gender vs app2_usage")
plot(app1$gender,app1$app_3,main="gender vs app3_usage")
plot(app1$gender,app1$app_4,main="gender vs app4_usage")
plot(app1$gender,app1$app_5,main="gender vs app5_usage")
plot(app1$gender,app1$app_6,main="gender vs app6_usage")
plot(app1$gender,app1$app_7,main="gender vs app7_usage")
plot(app1$gender,app1$app_8,main="gender vs app8_usage")
plot(app1$gender,app1$app_9,main="gender vs app9_usage")
plot(app1$gender,app1$app_10,main="gender vs app10_usage")

female<- app1[app1$gender == "F",]
male<- app1[app1$gender == "M",]
summary(female)
summary(male)
varfemale<- c(sd(female$app_1),sd(female$app_2),sd(female$app_3),sd(female$app_4),sd(female$app_5),sd(female$app_6),sd(female$app_4),sd(female$app_8),sd(female$app_9),sd(female$app_10))
varfemale
varmale<- c(sd(male$app_1),sd(male$app_2),sd(male$app_3),sd(male$app_4),sd(male$app_5),sd(male$app_6),sd(male$app_4),sd(male$app_8),sd(male$app_9),sd(male$app_10))
varmale

#covariance of apps
cor(app1[,3:12])

#(b)
#generate the logistic regression
M<- gsub("M","1",app1$gender)
gender<- gsub("F","0",M)
gender<-as.numeric(gender)
app.glm<- glm(gender ~ app1$app_1 + app1$app_2 + app1$app_3+ app1$app_4+ app1$app_5+ app1$app_6+ app1$app_7+ app1$app_8+ app1$app_9+ app1$app_10, family = binomial)
summary(app.glm)
app.glm1<- glm(gender ~  app1$app_1+ app1$app_2+app1$app_4+ app1$app_5+ app1$app_6+  app1$app_10, family = binomial)
colinear<- lm(app_1~app_2+app_4+ app_5+ app_6+app_10,data=app1)
summary(colinear)

#stepwise
ackwards = step(app.glm)

#(c)
# find critical value to classify female and male.
final.glm<- app.glm1
a1<- sort(unique(app1$app_1))
a2<- sort(unique(app1$app_2))
a4<- sort(unique(app1$app_4))
a5<- sort(unique(app1$app_5))
a6<- sort(unique(app1$app_6))
a10<- sort(unique(app1$app_10))
pred<- data.frame(a1,a2,a4,a5,a6,a10)
phat<- predict.glm(final.glm,pred,type="response")
oldapp <- cbind(phat,app1)
par(mfrow=c(1,1))
plot(oldapp$gender,oldapp$phat,main="TRUE Prediction")
male.pred<- oldapp[oldapp$gender == "M",]
female.pred<- oldapp[oldapp$gender == "F",]
t.test(male.pred$phat,female.pred$phat, paired=FALSE)

# data y simulated from a logistic regression model 
# with with three predictors, n=700
set.seed(750)
x = matrix(rnorm(4200),700,6)
lp = 1.05E-02 * x[,1] + 6.46E-02*x[,2] +1.57E-02*x[,3] -5.98E-03*x[,4] -2.96E-02*x[,5] -5.47E-02*x[,6]
p = 1/(1+exp(-lp))
y = runif(700)<p

# fit a logistic regression model
mod = glm(y~x[,1]+x[,2]+x[,3]+x[,4]+x[,5]+x[,6],family="binomial")

# using a cutoff of cut, calculate sensitivity, specificity, and classification rate
perf = function(cut, mod, y)
{
  yhat = (mod$fit>cut)
  w = which(y==1)
  sensitivity = mean( yhat[w] == 1 ) 
  specificity = mean( yhat[-w] == 0 ) 
  c.rate = mean( y==yhat ) 
  d = cbind(sensitivity,specificity)-c(1,1)
  d = sqrt( d[1]^2 + d[2]^2 ) 
  out = t(as.matrix(c(sensitivity, specificity, c.rate,d)))
  colnames(out) = c("sensitivity", "specificity", "c.rate", "distance")
  return(out)
}

s = seq(.01,.99,length=1000)
OUT = matrix(0,1000,4)
for(i in 1:1000) OUT[i,]=perf(s[i],mod,y)
par(mfrow=c(1,1))
plot(s,OUT[,1],xlab="Cutoff",ylab="Value",cex.lab=1.5,cex.axis=1.5,ylim=c(0,1),type="l",lwd=2,axes=FALSE,col=2)
axis(1,seq(0,1,length=5),seq(0,1,length=5),cex.lab=1.5)
axis(2,seq(0,1,length=5),seq(0,1,length=5),cex.lab=1.5)
lines(s,OUT[,2],col="darkgreen",lwd=2)
lines(s,OUT[,3],col=4,lwd=2)
lines(s,OUT[,4],col="darkred",lwd=2)
box()
legend(0,.35,col=c(2,"darkgreen",4,"darkred"),lwd=c(2,2,2,2),c("Sensitivity","Specificity","Classification Rate","Distance"))
abline(v=0.5)

#prediction
app2<-read.csv("MAT750.app2.csv")
final.glm<- app.glm1
a1<- sort(unique(app2$app_1))
a2<- sort(unique(app2$app_2))
a4<- sort(unique(app2$app_4))
a5<- sort(unique(app2$app_5))
a6<- sort(unique(app2$app_6))
a10<- sort(unique(app2$app_10))
pred<- data.frame(a1 =rep(a1,7),a2 =rep(a2,7),a4 =rep(a4,7),a5 =rep(a5,7),a6 =rep(a6,7),a10 =rep(a10,7))
pred<- pred[c(1:700),]
phat<- predict.glm(final.glm,pred,type="response")
phat<-phat[c(1:300)]
newapp <- cbind(phat,app2)
#classify gender
gender<-rep(0,300)
for (i in 1:300){
  if (newapp$phat[i] > 0.5)
    gender[i] <- 1
  else
    gender[i] <- 0
}
M<- gsub("1","M",gender)
gender.new<- gsub("0","F",M)
new.pred<- cbind(gender.new,app2)
summary(new.pred)

#(d)
#newdata summary
par(mfrow=c(2,5))
plot(new.pred$gender,new.pred$app_1,main="gender vs app1_usage")
plot(new.pred$gender,new.pred$app_2,main="gender vs app2_usage")
plot(new.pred$gender,new.pred$app_3,main="gender vs app3_usage")
plot(new.pred$gender,new.pred$app_4,main="gender vs app4_usage")
plot(new.pred$gender,new.pred$app_5,main="gender vs app5_usage")
plot(new.pred$gender,new.pred$app_6,main="gender vs app6_usage")
plot(new.pred$gender,new.pred$app_7,main="gender vs app7_usage")
plot(new.pred$gender,new.pred$app_8,main="gender vs app8_usage")
plot(new.pred$gender,new.pred$app_9,main="gender vs app9_usage")
plot(new.pred$gender,new.pred$app_10,main="gender vs app10_usage")

female.new<- new.pred[new.pred$gender == "F",]
male.new<- new.pred[new.pred$gender == "M",]
summary(female.new)
summary(male.new)
varfemale.new<- c(sd(female.new$app_1),sd(female.new$app_2),sd(female.new$app_3),sd(female.new$app_4),sd(female.new$app_5),sd(female.new$app_6),sd(female.new$app_4),sd(female.new$app_8),sd(female.new$app_9),sd(female.new$app_10))
varfemale.new
varmale.new<- c(sd(male.new$app_1),sd(male.new$app_2),sd(male.new$app_3),sd(male.new$app_4),sd(male.new$app_5),sd(male.new$app_6),sd(male.new$app_4),sd(male.new$app_8),sd(male.new$app_9),sd(male.new$app_10))
varmale.new

#3.(b)
#u1= mean(month)
#u2= mean(week)
#delta= Margin of non-inferiority
#kappa= samplesize(month) / sample size (week) = 1
#library(MVR)
#sd= pool.sd(month,week)
#alpha=0.05
#beta=0.20
#nweek = nmonth =(1+1/kappa)*(sd*(qnorm(1-alpha)+qnorm(1-beta))/(u1-u2-delta))^2)
#ceiling(nweek) 
