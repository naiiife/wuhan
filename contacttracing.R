
data <- read.table('contacttracing.csv',sep=',',head=T)
N = nrow(data)
data$ID <- 1:N
index = !is.na(data$是否感染他人.GRTR)
data1 <- data[index,]
index = data1$是否感染他人.GRTR==1 & data1$确诊时所在省份.SZSF!='海外'
data2 <- data1[index,]
N2 = nrow(data2)
infectionlist1 <- NULL
infectionlist2 <- NULL
infectionlist3 <- NULL
infectionlist4 <- NULL
infectionlist5 <- NULL
infectionlist6 <- NULL
infectionlist7 <- NULL
infectionlist8 <- NULL
infectionlist9 <- NULL

for (i in 1:N2){
if (!is.na(data2$被感染人病例编号..省份官方编号.1.GSBH1[i])){
chaini = data2[i,c('ID','确诊时所在省份.SZSF','确诊省份官方编号.GFBH',
        '被感染人病例编号..省份官方编号.1.GSBH1','感染日期.开始.1.RRQIDT1',
        '感染日期.结束.1.RRQDT1','首次出现症状日期.CXZZDT')]
infectionlist1 <- rbind(infectionlist1, chaini)
}
}
for (i in 1:N2){
if (!is.na(data2$被感染人病例编号..省份官方编号.2.GSBH2[i])){
chaini = data2[i,c('ID','确诊时所在省份.SZSF','确诊省份官方编号.GFBH',
        '被感染人病例编号..省份官方编号.2.GSBH2','感染日期.开始.2.RRQIDT2',
        '感染日期.结束.2.RRQDT2','首次出现症状日期.CXZZDT')]
infectionlist2 <- rbind(infectionlist2, chaini)
}
}
for (i in 1:N2){
if (!is.na(data2$被感染人病例编号..省份官方编号.3.GSBH3[i])){
chaini = data2[i,c('ID','确诊时所在省份.SZSF','确诊省份官方编号.GFBH',
        '被感染人病例编号..省份官方编号.3.GSBH3','感染日期.开始.3.RRQIDT3',
        '感染日期.结束.3.RRQDT3','首次出现症状日期.CXZZDT')]
infectionlist3 <- rbind(infectionlist3, chaini)
}
}
for (i in 1:N2){
if (!is.na(data2$被感染人病例编号..省份官方编号.4.GSBH4[i])){
chaini = data2[i,c('ID','确诊时所在省份.SZSF','确诊省份官方编号.GFBH',
        '被感染人病例编号..省份官方编号.4.GSBH4','感染日期.开始.4.RRQIDT4',
        '感染日期.结束.4.RRQDT4','首次出现症状日期.CXZZDT')]
infectionlist4 <- rbind(infectionlist4, chaini)
}
}
for (i in 1:N2){
if (!is.na(data2$被感染人病例编号..省份官方编号.5.GSBH5[i])){
chaini = data2[i,c('ID','确诊时所在省份.SZSF','确诊省份官方编号.GFBH',
        '被感染人病例编号..省份官方编号.5.GSBH5','感染日期.开始.5.RRQIDT5',
        '感染日期.结束.5.RRQDT5','首次出现症状日期.CXZZDT')]
infectionlist5 <- rbind(infectionlist5, chaini)
}
}
for (i in 1:N2){
if (!is.na(data2$被感染人病例编号..省份官方编号.6.GSBH6[i])){
chaini = data2[i,c('ID','确诊时所在省份.SZSF','确诊省份官方编号.GFBH',
        '被感染人病例编号..省份官方编号.6.GSBH6','感染日期.开始.6.RRQIDT6',
        '感染日期.结束.6.RRQDT6','首次出现症状日期.CXZZDT')]
infectionlist6 <- rbind(infectionlist6, chaini)
}
}
for (i in 1:N2){
if (!is.na(data2$被感染人病例编号..省份官方编号.7.GSBH7[i])){
chaini = data2[i,c('ID','确诊时所在省份.SZSF','确诊省份官方编号.GFBH',
        '被感染人病例编号..省份官方编号.7.GSBH7','感染日期.开始.7.RRQIDT7',
        '感染日期.结束.7.RRQDT7','首次出现症状日期.CXZZDT')]
infectionlist7 <- rbind(infectionlist7, chaini)
}
}
for (i in 1:N2){
if (!is.na(data2$被感染人病例编号..省份官方编号.8.GSBH8[i])){
chaini = data2[i,c('ID','确诊时所在省份.SZSF','确诊省份官方编号.GFBH',
        '被感染人病例编号..省份官方编号.8.GSBH8','感染日期.开始.8.RRQIDT8',
        '感染日期.结束.8.RRQDT8','首次出现症状日期.CXZZDT')]
infectionlist8 <- rbind(infectionlist8, chaini)
}
}
for (i in 1:N2){
if (!is.na(data2$被感染人病例编号..省份官方编号.9.GSBH9[i])){
chaini = data2[i,c('ID','确诊时所在省份.SZSF','确诊省份官方编号.GFBH',
        '被感染人病例编号..省份官方编号.9.GSBH9','感染日期.开始.9.RRQIDT9',
        '感染日期.结束.9.RRQDT9','首次出现症状日期.CXZZDT')]
infectionlist9 <- rbind(infectionlist9, chaini)
}
}

colnames(infectionlist1) <- c('ID','Province','InfectorID','InfecteeID',
     'InfectionB','InfectionE','Onsetor')
colnames(infectionlist2) <- c('ID','Province','InfectorID','InfecteeID',
     'InfectionB','InfectionE','Onsetor')
colnames(infectionlist3) <- c('ID','Province','InfectorID','InfecteeID',
     'InfectionB','InfectionE','Onsetor')
colnames(infectionlist4) <- c('ID','Province','InfectorID','InfecteeID',
     'InfectionB','InfectionE','Onsetor')
colnames(infectionlist5) <- c('ID','Province','InfectorID','InfecteeID',
     'InfectionB','InfectionE','Onsetor')
colnames(infectionlist6) <- c('ID','Province','InfectorID','InfecteeID',
     'InfectionB','InfectionE','Onsetor')
colnames(infectionlist7) <- c('ID','Province','InfectorID','InfecteeID',
     'InfectionB','InfectionE','Onsetor')
colnames(infectionlist8) <- c('ID','Province','InfectorID','InfecteeID',
     'InfectionB','InfectionE','Onsetor')
colnames(infectionlist9) <- c('ID','Province','InfectorID','InfecteeID',
     'InfectionB','InfectionE','Onsetor')

infectionlist <- rbind(infectionlist1,infectionlist2,infectionlist3,
infectionlist4,infectionlist5,infectionlist6,infectionlist7,
infectionlist8,infectionlist9)
N0 = nrow(infectionlist)
infectionlist <- cbind(infectionlist,rep(NA,N0),rep(NA,N0))
colnames(infectionlist) <- c('LineNOor','Province','InfectorID','InfecteeID',
     'InfectionB','InfectionE','Onsetor','Onsetee','LineNOee')

overing = rep(0,N0)
for (i in 1:N0){
index = which(data$确诊时所在省份.SZSF==infectionlist[i,2] & 
      data$确诊省份官方编号.GFBH==infectionlist[i,4])
if (length(index)==1) {
if (!is.na(data$是否有湖北接触史.SHJC[index])){
if (na.omit(data$是否有湖北接触史.SHJC[index])==1) overing[i]=1
}
infectionlist[i,9] = data$ID[index]
infectionlist[i,8] = as.character(data$首次出现症状日期.CXZZDT[index])
}
}
infectionlist <- infectionlist[overing==0,]
N0 = nrow(infectionlist)

overing = rep(0,N0)
for (i in 1:N0){
if (length(which(infectionlist[,9]==infectionlist[i,9]))>1) overing[i]=1
}
infectionlist <- infectionlist[overing==0,]
N0 = nrow(infectionlist)


Serial <- rep(NA,N0)
for (i in 1:N0){
onsetor = as.character(infectionlist[i,7])
onsetee = as.character(infectionlist[i,8])
if (onsetor!='' & onsetee!=''){
Serial[i] <- as.numeric(as.Date(onsetee)-as.Date(onsetor))
}
}
#write.table(infectionlist,'infectionlist.csv',sep=',')

serial <- na.omit(Serial)
hist(serial,breaks=20,main='Histogram of Serial Intervals',xlab='Days')
#savePlot('serial',type='png')
length(serial)
mean(serial)
sd(serial)
