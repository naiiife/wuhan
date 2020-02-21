#setwd("")
case.info <- read.csv(file = "contact-tracing.csv",header = T)
dim(case.info)
names(case.info)
tail(case.info)
case.info <- case.info[1:221,1:10]
str(case.info)
tail(case.info)
case.info$date.symptom <- as.Date(case.info$date.symptom,"%m/%d/%Y")

dat1 <- case.info[!is.na(case.info$id.infected),]
dat1$d2 <- "1/1/0000"
dat1$d2 <- as.Date(dat1$d2,"%m/%d/%Y")

index <- which(!is.na(case.info$id.infected))
id  <- case.info$id.infected[index]
pro <- case.info$province[index]

for (j in 1:length(index)){
  print(j)
  d <- case.info[case.info$ID==id[j]&case.info$province==pro[j],]$date.symptom
  if(!is.na(d)){dat1$d2[j] <- d}else{dat1$d2[j]=NA}
}

dat1$d2
dat1

interval <- dat1$date.symptom-dat1$d2
dat1[which(interval==0),]
dat1[which(interval<0),]

interval <- interval[!is.na(interval)]
interval <- as.numeric(interval)
interval <- as.integer(interval)
length(interval)
summary(interval)
table(interval)
hist(interval, main = "Histogram of Serial Interval",xlab = "Serial Interval",xlim=c(-6,14))


