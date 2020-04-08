## This software is licenced under the "Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International Public License" (http://creativecommons.org/licenses/by/4.0/)
## TAKE CHROMOPAINTER OUTPUT AND FORMS THE HAPLOTYPE SHARING PATTERNS OF A TARGET GROUP AS A MIXTURE OF THAT FROM A SET OF SURROGATE GROUPS, INFERRING THE MIXTURE COEFFICIENTS USING MARKOV-CHAIN-MONTE-CARLO (MCMC) 
## AS DESCRIBED IN "Latin Americans show wide-spread Converso ancestry and the imprint of local Native ancestry on physical appearance" (2018) by Chacon-Duque et al (preprint: https://www.biorxiv.org/content/early/2018/01/23/252155)
## HAS SLIGHTLY DIFFERENT PROPOSAL DENSITY TO THAT DESCRIBED IN THE PAPER LISTED ABOVE -- SEE INSTRUCTION MANUAL (THANKS TO IDA MOLTKE AND RYAN WAPLES FOR SPOTTING ERRORS)

## usage:   (1) R < sourcefindv2.R parameter.infile.list --no-save > output.out
 
## example:  R < sourcefindv2.R example/BrahuiYorubaSimulation.SourcefindParamfile.txt --no-save > output.out

########################################
## COMMAND LINE INPUT:

usage=function()
{
	print(noquote("run using: R < sourcefindv2.R parameter_infile --no-save > screen_output"))
	print(noquote("parameter input file format (NO defaults, all fields must be entered, even if not used):"))
	print(noquote("self.copy.ind: [0,1]"))
	print(noquote("num.surrogates: [2,...,S = number of surrogates]"))
	print(noquote("exp.num.surrogates: [1,2,...,S]"))
	print(noquote("input.file.ids: [input.filename1]"))
	print(noquote("input.file.copyvectors: [input.filename2]"))
	print(noquote("save.file: [output.filename]"))
	print(noquote("copyvector.popnames: [pop_1 pop_2 ... pop_j]"))
	print(noquote("surrogate.popnames: [pop_1 pop_2 ... pop_k]"))
	print(noquote("target.popnames: [pop_target_1 ... pop_target_m]"))
	print(noquote("num.slots: [1,2,...]"))
	print(noquote("num.iterations: [1,2,...]"))
	print(noquote("num.burnin: [1,2,...]"))
	print(noquote("num.thin: [1,2,...]"))
}

temp=commandArgs()

param.infile=as.character(temp[2])
if (param.infile=="help"){usage();q(save='no')}

######################################
## COMMAND LINE CHECK AND READ IN:

error.input.message=function(file.name)
  {
    print(paste("Something wrong with input file ",file.name,". See below. Exiting...",sep=''))
    usage()
    q(save='no')
  }
line.check=function(file.name,skip.val,match.val)
  {
    if (as.character(read.table(file.name,skip=skip.val,nrows=1,as.is=TRUE)[1]) != as.character(match.val))
      {
        error.input.message(file.name)
      }
    if (as.character(read.table(file.name,skip=skip.val,nrows=1,as.is=TRUE)[1]) == as.character(match.val))
      return(read.table(file.name,skip=skip.val,nrows=1,as.is=TRUE)[2:length(read.table(file.name,skip=skip.val,nrows=1,as.is=TRUE))])
  }

          ## line read in and checks:
self.copy.ind=as.integer(line.check(param.infile,0,"self.copy.ind:"))
if (length(self.copy.ind)!=1 || (self.copy.ind!=0 && self.copy.ind!=1)) error.input.message(param.infile)
num.pops=as.integer(line.check(param.infile,1,"num.surrogates:"))
if (length(num.pops)!=1 || num.pops<=1) error.input.message(param.infile)
surr.exp=as.integer(line.check(param.infile,2,"exp.num.surrogates:"))
if (length(surr.exp)!=1 || surr.exp<=0) error.input.message(param.infile)
id.file=as.character(line.check(param.infile,3,"input.file.ids:"))
if (length(id.file)!=1) error.input.message(param.infile)
copyvector.file=as.character(line.check(param.infile,4,"input.file.copyvectors:"))
if (length(copyvector.file)!=1) error.input.message(param.infile)
save.file.props=as.character(line.check(param.infile,5,"save.file:"))
if (length(save.file.props)!=1) error.input.message(param.infile)
donor.pops.all=unique(as.character(line.check(param.infile,6,"copyvector.popnames:")))
if (length(donor.pops.all)<=1) error.input.message(param.infile)
newnames=as.character(line.check(param.infile,7,"surrogate.popnames:"))
if (length(unique(newnames))<=1) error.input.message(param.infile)
recipient.pops=as.character(line.check(param.infile,8,"target.popnames:"))
if (length(recipient.pops)<=0) error.input.message(param.infile)
num.slots = as.integer(line.check(param.infile,9,"num.slots:"))
if (length(num.slots)!=1 || num.slots <= 1) error.input.message(param.infile)
num.runs = as.integer(line.check(param.infile,10,"num.iterations:"))
if (length(num.runs)!=1 || num.runs <= 1) error.input.message(param.infile)
burn.in = as.integer(line.check(param.infile,11,"num.burnin:"))
if (length(burn.in)!=1 || burn.in <= 1) error.input.message(param.infile)
thin.val = as.integer(line.check(param.infile,12,"num.thin:"))
if (length(thin.val)!=1 || thin.val <= 1) error.input.message(param.infile)

                           ## standard checks:
if (burn.in>=num.runs) {print("You have specified more burn-in iterations than total iterations. Exiting...."); q(save='no')}
if ((burn.in+thin.val)>num.runs) {print("Your specified burn-in iterations, total iterations, and thinning value will result in no samples. Exiting...."); q(save='no')}
if (surr.exp>num.pops) {print("Your specified mean number of surrogates per iteration is greater than the total number of allowed surrogates per iteration. Exiting...."); q(save='no')}

#############################################################
##############################################################
## (0) FUNCTIONS:
##############################################################
#############################################################

options(scipen=999)

multi.prob.func=function(x,target.vec) return(sum(target.vec*log(x)))

sourcefind=function(target.vec,surr.mat,num.pops,surr.exp,num.slots,num.runs,burn.in,thin.val,target.name,surr.pops,save.file)
{
	add.prior='yes'
	num.surr=dim(surr.mat)[1]          ## total number of surrogates in data
	if (num.pops>num.surr) num.pops=num.surr
	if (surr.exp>num.surr) surr.exp=num.surr
	if (num.slots<num.pops) num.pops=num.slots
	if (num.slots<surr.exp) surr.exp=num.slots
	pops.toreplace=1

        prior.prob.vec=(surr.exp^(1:num.pops)/factorial(1:num.pops))/sum(surr.exp^(1:num.pops)/factorial(1:num.pops))

        surr.tochoose.orig=as.integer(sample(as.character(1:num.surr),num.pops,replace=T))
	surr.tochoose=as.integer(sample(as.character(surr.tochoose.orig),num.slots,replace=T))
	prev.multinom=multi.prob.func(apply((1.0/num.slots)*matrix(surr.mat[surr.tochoose,],ncol=dim(surr.mat)[2]),2,sum),target.vec)
	if (add.prior=='yes') prev.multinom=prev.multinom+log(prior.prob.vec[length(unique(surr.tochoose))])
	surr.tochoose.prev=surr.tochoose
	surr.tochoose.final=multi.prob.final=NULL
	num.sample=0
	for (m in 1:num.runs)
	{
		#ind.to.replace=as.integer(sample(as.character(unique(surr.tochoose)),pops.toreplace,replace=T))
		ind.to.replace=as.integer(sample(as.character(sort(unique(surr.tochoose))),pops.toreplace,prob=table(sort(unique(surr.tochoose)))/num.slots,replace=T))
		for (i in 1:pops.toreplace)
		{
			num.slots.toreplace.final=sample(1:sum(surr.tochoose==ind.to.replace[i]),1)
			random.val1=runif(1)
			if (random.val1<0.1) num.slots.toreplace.final=sum(surr.tochoose==ind.to.replace[i])
			slots.to.replace=as.integer(sample(as.character((1:num.slots)[surr.tochoose==ind.to.replace[i]]),num.slots.toreplace.final,replace=F))
			slots.to.keep=(1:num.slots)[-slots.to.replace]
			all.surr=integer(0)
			if (length(slots.to.keep)>0) all.surr=unique(surr.tochoose[slots.to.keep])
			if (length(all.surr)<num.pops && length((1:num.surr)[-all.surr])>0) all.surr=c(all.surr,as.integer(sample(as.character((1:num.surr)[-all.surr]),num.pops-length(all.surr),replace=F)))
			other.surr=setdiff(all.surr,ind.to.replace[i])
			if (random.val1<0.1) other.surr=setdiff(1:num.surr,c(unique(surr.tochoose[slots.to.keep]),ind.to.replace[i]))
			if (length(other.surr)==0)
			{
				other.surr=as.integer(sample(as.character((1:num.surr)[-ind.to.replace[i]]),1))
				slots.to.replace=1:num.slots
			}
			surr.tochoose[slots.to.replace]=rep(as.integer(sample(as.character(other.surr),1)),length(slots.to.replace))
 		}
			  ## calculate new prob, and accept if better:
		new.multinom=multi.prob.func(apply((1.0/num.slots)*matrix(surr.mat[surr.tochoose,],ncol=dim(surr.mat)[2]),2,sum),target.vec)
		if (add.prior=='yes') new.multinom=new.multinom+log(prior.prob.vec[length(unique(surr.tochoose))])
		accept.prob.multinom=min(c(1,exp(new.multinom-prev.multinom)))
		random.val=runif(1)
		if (random.val<=accept.prob.multinom)
		{
			prev.multinom=new.multinom
			surr.tochoose.prev=surr.tochoose
		}
		if (random.val>accept.prob.multinom)
		{
			new.multinom=prev.multinom
			surr.tochoose=surr.tochoose.prev
		}
		if (m==(burn.in+1+num.sample*thin.val))
		{
			multi.prob.final=c(multi.prob.final,new.multinom)
			table.surr=table(surr.tochoose)
			surr.tochoose.new=rep(0,num.surr)
			surr.tochoose.new[as.integer(names(table.surr))]=table.surr/num.slots
			surr.tochoose.final=rbind(surr.tochoose.final,surr.tochoose.new)
			num.sample=num.sample+1
		}
     }
     surr.tochoose.toprint=matrix(0,nrow=dim(surr.tochoose.final)[1],ncol=length(surr.pops))
     surr.tochoose.toprint[,match(rownames(surr.mat),surr.pops)]=surr.tochoose.final
     ###return(surr.tochoose.final) 
     ###return(surr.tochoose.final[multi.prob.final==max(multi.prob.final),][1,])
     write.table(rbind(c("target","posterior.prob",surr.pops),cbind(target.name,multi.prob.final,surr.tochoose.toprint)),file=save.file,row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
     return(1)
}


##############################################################
#############################################################
### (I) READ IN DATA AND MAKE COPY-VECTORS:
##############################################################
#############################################################

donor.pops.all2=unique(newnames)

                        ## GET IDS FILE:
id.mat=read.table(id.file,as.is=TRUE)
to.remove.id=(1:dim(id.mat)[1])[as.character(id.mat[,3])=='0']
if (length(to.remove.id)>0) id.mat=id.mat[-to.remove.id,]
if (length(dim(id.mat))==0)
{
	print(paste("SOMETHING WRONG WITH ",id.file," -- NO NON-EXCLUDED INDS? Exiting....",sep=''))
	q(save='no')
}

                      ## GET RAW COPY PROPS FOR ALL DONORS:
probs=read.table(copyvector.file,header=TRUE,check.names=FALSE)
rownames.copyvector=as.character(probs[,1])
colnames.copyvector=as.character(read.table(copyvector.file,as.is=TRUE,nrows=1,check.names=FALSE)[-1])
for (i in 1:length(donor.pops.all))
{
	if (!is.element(donor.pops.all[i],as.character(id.mat[,2])) && !is.element(donor.pops.all[i],colnames.copyvector))
	{
		print(paste("COPY VECTOR COLUMN LABEL ",donor.pops.all[i]," NOT FOUND IN ",id.file," OR COLUMNS OF ",copyvector.file,"! Exiting....",sep=''))
		q(save='no')
	}
}
for (i in 1:length(donor.pops.all2))
{
	if (!is.element(donor.pops.all2[i],as.character(id.mat[,2])) && !is.element(donor.pops.all2[i],rownames.copyvector))
	{
		print(paste("SURROGATE POPULATION ",donor.pops.all2[i]," NOT FOUND IN ",id.file," OR ROWS OF ",copyvector.file,"! Exiting....",sep=''))
		q(save='no')
	}
}
for (i in 1:length(recipient.pops))
{
	if (!is.element(recipient.pops[i],as.character(id.mat[,2])) && !is.element(recipient.pops[i],rownames.copyvector))
	{
		print(paste("TARGET POPULATION ",recipient.pops[i]," NOT FOUND IN ",id.file," OR ROWS OF ",copyvector.file,"! Exiting....",sep=''))
		q(save='no')
	}
}

                                   ## combine columns across copy-vector pops:
predmat.orig=NULL
for (i in 1:length(donor.pops.all))
  {
    id.labels.i=c(donor.pops.all[i],as.character(id.mat[,1])[as.character(id.mat[,2])==donor.pops.all[i]])
    match.i=NULL
    for (j in 1:length(id.labels.i)) match.i=c(match.i,which(as.character(colnames.copyvector)==id.labels.i[j]))
    if (length(match.i)==0)
    {
	print(paste("NO INDS OF ",donor.pops.all[i]," FOUND AMONG COLUMNS OF ",copyvector.file,"! Exiting....",sep=''))
	q(save='no')
    }
    predmat.orig=cbind(predmat.orig,apply(matrix(as.matrix(probs[,2:dim(probs)[2]][,match.i]),nrow=dim(probs)[1]),1,sum))
  }
rownames(predmat.orig)=rownames.copyvector
colnames(predmat.orig)=donor.pops.all

                                   ## combine rows across donor pops:
predmat=NULL
for (i in 1:length(donor.pops.all2))
  {
    id.labels.i=c(donor.pops.all2[i],as.character(id.mat[,1])[as.character(id.mat[,2])==donor.pops.all2[i]])
    match.i=NULL
    for (j in 1:length(id.labels.i)) match.i=c(match.i,which(as.character(rownames.copyvector)==id.labels.i[j]))
    if (length(match.i)==0)
    {
	print(paste("NO INDS OF ",donor.pops.all2[i]," FOUND AMONG ROWS OF ",copyvector.file,"! Exiting....",sep=''))
	q(save='no')
    }
    predmat=rbind(predmat,apply(matrix(as.matrix(predmat.orig[match.i,]),ncol=dim(predmat.orig)[2]),2,mean))
  }
rownames(predmat)=donor.pops.all2
colnames(predmat)=donor.pops.all

                       ## GET RAW COPY PROPS FOR RECIPIENT:
recipient.mat=NULL
for (i in 1:length(recipient.pops))
  {
    id.labels.i=c(recipient.pops[i],as.character(id.mat[,1])[as.character(id.mat[,2])==recipient.pops[i]])
    match.i=NULL
    for (j in 1:length(id.labels.i)) match.i=c(match.i,which(as.character(rownames.copyvector)==id.labels.i[j]))
    if (length(match.i)==0)
    {
	print(paste("NO INDS OF ",recipient.pops[i]," FOUND AMONG ROWS OF ",copyvector.file,"! Exiting....",sep=''))
	q(save='no')
    }
    recipient.mat=rbind(recipient.mat,apply(matrix(as.matrix(predmat.orig[match.i,]),ncol=dim(predmat.orig)[2]),2,mean))
  }
rownames(recipient.mat)=recipient.pops
colnames(recipient.mat)=donor.pops.all

                        ## run SOURCEFIND on each recipient pop:
write("###INFERRED PROPORTIONS OF ANCESTRY",file=save.file.props,ncolumns=1)
for (i in 1:length(recipient.pops))
{
	print(paste("Analysing target ",i," of ",length(recipient.pops)," -- ",recipient.pops[i],".....",sep=''))
	target.vec=recipient.mat[i,]
	surr.mat=predmat
	surr.pops=rownames(predmat)
	row.i=(1:length(donor.pops.all2))[donor.pops.all2==recipient.pops[i]]
	col.i=(1:length(donor.pops.all))[donor.pops.all==recipient.pops[i]]
	if ((self.copy.ind==0 && length(row.i)>0) || (self.copy.ind==1 && length(col.i)==0 && length(row.i)>0)) surr.mat=surr.mat[-row.i,]
	if (self.copy.ind==1 && length(col.i)>0)
	{
		self.copy.vec=rep(0,length(donor.pops.all))
		self.copy.vec[col.i]=1
		if (length(row.i)>0) surr.mat[row.i,]=self.copy.vec
		if (length(row.i)==0)
		{
			surr.mat=rbind(surr.mat,self.copy.vec)
			surr.pops=c(surr.pops,recipient.pops[i])
			rownames(surr.mat)=surr.pops
		}
	}
	for (j in 1:dim(surr.mat)[1]) surr.mat[j,]=surr.mat[j,]/sum(surr.mat[j,])
	sourcefind.run=sourcefind(target.vec,surr.mat,num.pops,surr.exp,num.slots,num.runs,burn.in,thin.val,recipient.pops[i],surr.pops,save.file.props)
}

print("Finished!")

q(save='no')
