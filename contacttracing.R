
data <- read.table('contacttracing.csv',sep=',',head=T)
N = nrow(data)
data$ID <- 1:N
index = !is.na(data$�Ƿ��Ⱦ����.GRTR)
data1 <- data[index,]
index = data1$�Ƿ��Ⱦ����.GRTR==1 & data1$ȷ��ʱ����ʡ��.SZSF!='����'
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
if (!is.na(data2$����Ⱦ�˲������..ʡ�ݹٷ����.1.GSBH1[i])){
chaini = data2[i,c('ID','ȷ��ʱ����ʡ��.SZSF','ȷ��ʡ�ݹٷ����.GFBH',
        '����Ⱦ�˲������..ʡ�ݹٷ����.1.GSBH1','��Ⱦ����.��ʼ.1.RRQIDT1',
        '��Ⱦ����.����.1.RRQDT1','�״γ���֢״����.CXZZDT')]
infectionlist1 <- rbind(infectionlist1, chaini)
}
}
for (i in 1:N2){
if (!is.na(data2$����Ⱦ�˲������..ʡ�ݹٷ����.2.GSBH2[i])){
chaini = data2[i,c('ID','ȷ��ʱ����ʡ��.SZSF','ȷ��ʡ�ݹٷ����.GFBH',
        '����Ⱦ�˲������..ʡ�ݹٷ����.2.GSBH2','��Ⱦ����.��ʼ.2.RRQIDT2',
        '��Ⱦ����.����.2.RRQDT2','�״γ���֢״����.CXZZDT')]
infectionlist2 <- rbind(infectionlist2, chaini)
}
}
for (i in 1:N2){
if (!is.na(data2$����Ⱦ�˲������..ʡ�ݹٷ����.3.GSBH3[i])){
chaini = data2[i,c('ID','ȷ��ʱ����ʡ��.SZSF','ȷ��ʡ�ݹٷ����.GFBH',
        '����Ⱦ�˲������..ʡ�ݹٷ����.3.GSBH3','��Ⱦ����.��ʼ.3.RRQIDT3',
        '��Ⱦ����.����.3.RRQDT3','�״γ���֢״����.CXZZDT')]
infectionlist3 <- rbind(infectionlist3, chaini)
}
}
for (i in 1:N2){
if (!is.na(data2$����Ⱦ�˲������..ʡ�ݹٷ����.4.GSBH4[i])){
chaini = data2[i,c('ID','ȷ��ʱ����ʡ��.SZSF','ȷ��ʡ�ݹٷ����.GFBH',
        '����Ⱦ�˲������..ʡ�ݹٷ����.4.GSBH4','��Ⱦ����.��ʼ.4.RRQIDT4',
        '��Ⱦ����.����.4.RRQDT4','�״γ���֢״����.CXZZDT')]
infectionlist4 <- rbind(infectionlist4, chaini)
}
}
for (i in 1:N2){
if (!is.na(data2$����Ⱦ�˲������..ʡ�ݹٷ����.5.GSBH5[i])){
chaini = data2[i,c('ID','ȷ��ʱ����ʡ��.SZSF','ȷ��ʡ�ݹٷ����.GFBH',
        '����Ⱦ�˲������..ʡ�ݹٷ����.5.GSBH5','��Ⱦ����.��ʼ.5.RRQIDT5',
        '��Ⱦ����.����.5.RRQDT5','�״γ���֢״����.CXZZDT')]
infectionlist5 <- rbind(infectionlist5, chaini)
}
}
for (i in 1:N2){
if (!is.na(data2$����Ⱦ�˲������..ʡ�ݹٷ����.6.GSBH6[i])){
chaini = data2[i,c('ID','ȷ��ʱ����ʡ��.SZSF','ȷ��ʡ�ݹٷ����.GFBH',
        '����Ⱦ�˲������..ʡ�ݹٷ����.6.GSBH6','��Ⱦ����.��ʼ.6.RRQIDT6',
        '��Ⱦ����.����.6.RRQDT6','�״γ���֢״����.CXZZDT')]
infectionlist6 <- rbind(infectionlist6, chaini)
}
}
for (i in 1:N2){
if (!is.na(data2$����Ⱦ�˲������..ʡ�ݹٷ����.7.GSBH7[i])){
chaini = data2[i,c('ID','ȷ��ʱ����ʡ��.SZSF','ȷ��ʡ�ݹٷ����.GFBH',
        '����Ⱦ�˲������..ʡ�ݹٷ����.7.GSBH7','��Ⱦ����.��ʼ.7.RRQIDT7',
        '��Ⱦ����.����.7.RRQDT7','�״γ���֢״����.CXZZDT')]
infectionlist7 <- rbind(infectionlist7, chaini)
}
}
for (i in 1:N2){
if (!is.na(data2$����Ⱦ�˲������..ʡ�ݹٷ����.8.GSBH8[i])){
chaini = data2[i,c('ID','ȷ��ʱ����ʡ��.SZSF','ȷ��ʡ�ݹٷ����.GFBH',
        '����Ⱦ�˲������..ʡ�ݹٷ����.8.GSBH8','��Ⱦ����.��ʼ.8.RRQIDT8',
        '��Ⱦ����.����.8.RRQDT8','�״γ���֢״����.CXZZDT')]
infectionlist8 <- rbind(infectionlist8, chaini)
}
}
for (i in 1:N2){
if (!is.na(data2$����Ⱦ�˲������..ʡ�ݹٷ����.9.GSBH9[i])){
chaini = data2[i,c('ID','ȷ��ʱ����ʡ��.SZSF','ȷ��ʡ�ݹٷ����.GFBH',
        '����Ⱦ�˲������..ʡ�ݹٷ����.9.GSBH9','��Ⱦ����.��ʼ.9.RRQIDT9',
        '��Ⱦ����.����.9.RRQDT9','�״γ���֢״����.CXZZDT')]
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
index = which(data$ȷ��ʱ����ʡ��.SZSF==infectionlist[i,2] & 
      data$ȷ��ʡ�ݹٷ����.GFBH==infectionlist[i,4])
if (length(index)==1) {
if (!is.na(data$�Ƿ��к����Ӵ�ʷ.SHJC[index])){
if (na.omit(data$�Ƿ��к����Ӵ�ʷ.SHJC[index])==1) overing[i]=1
}
infectionlist[i,9] = data$ID[index]
infectionlist[i,8] = as.character(data$�״γ���֢״����.CXZZDT[index])
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