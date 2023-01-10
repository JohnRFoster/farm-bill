############################################################################
###
### Data processing functions for multi-method removal data
###
### January 7, 2022
### Amy J. Davis
###
############################################################################






###
### Function to calculate the area impacted by removal efforts at each time and removal event
###

area.of.impact.finder<-function(rdata,gtypes,ymat,ardf,GPS=TRUE,detdf,plotareas=FALSE){
  require(sp)
  require(rgeos)
  rdata$GroupMth=paste(rdata$Site,rdata$MYNum,sep=".")
  rdata$GroupMthDay=paste(rdata$Site,rdata$MYNum,rdata$Day,sep="-")
  rdata$nm=as.numeric(factor(rdata$GroupMth,ordered=TRUE,levels=unique(rdata$GroupMth)))
  rdata$MthYear=format(rdata$Date,"%b-%Y")

  agb=matrix(detdf[match(gtypes,detdf$Method),"Area"],length(unique(rdata$GroupMth)),dim(gtypes)[2])



  if(GPS==TRUE){

    ###
    ### Processing to get the proportion of the area impacted by removal pass
    ###
    ### Creating unique site, month, and event unifiers

    ### Getting the number of each event by group and month
    tfv=rdata[,c("GroupMthDay","GroupMth")]
    tfv=tfv[!duplicated(tfv),]
    tfv$nm=as.numeric(factor(tfv$GroupMth,ordered=TRUE,levels=unique(tfv$GroupMth)))
    tfv$Dayn=ave(tfv$nm,tfv$nm,FUN=seq_along)
    rdata$Dayn=tfv[match(rdata$GroupMthDay,tfv$GroupMthDay),"Dayn"]

    ardfgm=data.frame(SiteMth=0,Day=0,Dayn=0,nm=0,Area=0)

    ### Slow way to get an visualize area of impact by event
    for(k in unique(rdata$GroupMthDay)){
      kind=which(unique(rdata$GroupMthDay)==k)
      kdat=rdata[rdata$GroupMthDay==k,]
      spdf=SpatialPointsDataFrame(coords = rdata[rdata$GroupMthDay==k,c("Longitude","Latitude")],
                                  data = data.frame(1:dim(rdata[rdata$GroupMthDay==k,])[1]))
      sppoly=gConvexHull(spdf,byid = FALSE)
      spbuff=gBuffer(spgeom = sppoly,byid=FALSE,width=0.02)
      ardfgm[kind,"SiteMth"]=kdat$GroupMth[1]
      ardfgm[kind,"Day"]=kdat$Day[1]
      ardfgm[kind,"Dayn"]=kdat$Dayn[1]
      ardfgm[kind,"nm"]=kdat$nm[1]
      ardfgm[kind,"Area"]=gArea(spbuff)/0.00008

      if(plotareas==TRUE){
        par(mfrow=c(3,3))

        plot(rdata$Longitude[rdata$GroupMthDay==k],rdata$Latitude[rdata$GroupMthDay==k],pch=16,
             xlim=c(min(rdata$Longitude)*1.00019,max(rdata$Longitude)*0.9998),
             ylim=c(min(rdata$Latitude)*0.9992,max(rdata$Latitude)*1.0009),
             xlab="Longitude",ylab="Latitude")
        plot(spbuff,add=TRUE)
        title(paste("Group ",(rdata$MthYear)[which(unique(rdata$GroupMthDay)==k)]))
      }

    }
    agbe=agb
    agbe[cbind(ardfgm$nm,ardfgm$Dayn)]=ardfgm$Area
    agbe=agbe/agb

    gbea=gbe/agbe
    gbea[is.na(gbea)]=0

    Areai=as.matrix(agb)
    Areaper=matrix(pmin(1,Areai/rep(areatot,as.numeric(table(ymat[,1])))),dim(Areai)[1],dim(Areai)[2])
    Areaper[is.na(Areaper)]=0

    par(mfrow=c(1,1),mar=c(4,5,1,1)+0.1)
  }else if(GPS==FALSE){
    areatots=matrix(ardf[match(ymat$Site,ardf$Site),"Area"],length(unique(rdata$GroupMth)),dim(gtypes)[2],byrow = FALSE)
    Areaper=agb^(log(areatots*((gbe*agb)/(areatots+(gbe*agb))))/log(agb))/areatots
    Areaper[is.na(Areaper)]=0

  }
  Areaper

}




###################
###
### Editing the area finder function to work with primary and secondary period coding instead of month and day coding
###
###################


###
### Function to calculate the area impacted by removal efforts at each time and removal event
###

area.of.impact.finder2<-function(rdata,gtypes,ymat,ardf,gbe,areatot,GPS=TRUE,detdf,plotareas=FALSE){
  require(sp)
  require(rgeos)
  rdata$GroupMth=paste(rdata$Site,rdata$PPNum,sep=".")
  rdata$GroupMthDay=paste(rdata$Site,rdata$PPNum,rdata$SPNum,sep="-")
  rdata$nm=as.numeric(factor(rdata$GroupMth,ordered=TRUE,levels=unique(rdata$GroupMth)))
  rdata$MthYear=format(ymd_hms(rdata$Date),"%b-%Y")

  agb=matrix(detdf[match(gtypes,detdf$Method),"Area"],length(unique(rdata$GroupMth)),dim(gtypes)[2])



  if(GPS==TRUE){

    ###
    ### Processing to get the proportion of the area impacted by removal pass
    ###
    ### Creating unique site, month, and event unifiers

    ### Getting the number of each event by group and month
    tfv=rdata[,c("GroupMthDay","GroupMth")]
    tfv=tfv[!duplicated(tfv),]
    tfv$nm=as.numeric(factor(tfv$GroupMth,ordered=TRUE,levels=unique(tfv$GroupMth)))
    tfv$SPNumn=ave(tfv$nm,tfv$nm,FUN=seq_along)
    rdata$SPNumn=tfv[match(rdata$GroupMthDay,tfv$GroupMthDay),"SPNumn"]

    ardfgm=data.frame(SiteMth=0,SPNum=0,SPNumn=0,nm=0,Area=0)

    ### Slow way to get an visualize area of impact by event
    for(k in unique(rdata$GroupMthDay)){
      kind=which(unique(rdata$GroupMthDay)==k)
      kdat=rdata[rdata$GroupMthDay==k,]
      spdf=SpatialPointsDataFrame(coords = rdata[rdata$GroupMthDay==k,c("Longitude","Latitude")],
                                  data = data.frame(1:dim(rdata[rdata$GroupMthDay==k,])[1]))
      sppoly=gConvexHull(spdf,byid = FALSE)
      spbuff=gBuffer(spgeom = sppoly,byid=FALSE,width=0.02)
      ardfgm[kind,"SiteMth"]=kdat$GroupMth[1]
      ardfgm[kind,"SPNum"]=kdat$SPNum[1]
      ardfgm[kind,"SPNumn"]=kdat$SPNumn[1]
      ardfgm[kind,"nm"]=kdat$nm[1]
      ardfgm[kind,"Area"]=gArea(spbuff)/0.00008

      if(plotareas==TRUE){
        par(mfrow=c(3,3))

        plot(rdata$Longitude[rdata$GroupMthDay==k],rdata$Latitude[rdata$GroupMthDay==k],pch=16,
             xlim=c(min(rdata$Longitude)*1.00019,max(rdata$Longitude)*0.9998),
             ylim=c(min(rdata$Latitude)*0.9992,max(rdata$Latitude)*1.0009),
             xlab="Longitude",ylab="Latitude")
        plot(spbuff,add=TRUE)
        title(paste("Group ",(rdata$MthYear)[which(unique(rdata$GroupMthDay)==k)]))
      }

    }
    agbe=agb
    agbe[cbind(ardfgm$nm,ardfgm$SPNumn)]=ardfgm$Area
    agbe=agbe/agb

    gbea=gbe/agbe
    gbea[is.na(gbea)]=0

    Areai=as.matrix(agb)
    Areaper=matrix(pmin(1,Areai/rep(areatot,as.numeric(table(ymat[,1])))),dim(Areai)[1],dim(Areai)[2])
    Areaper[is.na(Areaper)]=0

    par(mfrow=c(1,1),mar=c(4,5,1,1)+0.1)
  }else if(GPS==FALSE){
    areatots=matrix(ardf[match(ymat$Site,ardf$Site),"Area"],length(unique(rdata$GroupMth)),dim(gtypes)[2],byrow = FALSE)
    Areaper=agb^(log(areatots*((gbe*agb)/(areatots+(gbe*agb))))/log(agb))/areatots
    Areaper[is.na(Areaper)]=0

  }
  Areaper

}

###################
###
### Creating a flexible primary period sequence based on a set number of months, default is 1 month
###
###################

primaryperiodsequence<-function(PPmonths=1,rdate=rdata$Date){
  ceiling(((as.numeric(format(rdate,"%m"))+(as.numeric(format(rdate,"%Y"))-min(as.numeric(format(rdate,"%Y"))))*12)-min(as.numeric(format(rdate,"%m"))+(as.numeric(format(rdate,"%Y"))-min(as.numeric(format(rdate,"%Y"))))*12)+1)/PPmonths)
}


###################
###
### Creating a flexible secondary period sequence based on a set number of days, default is 1 day
###
###################
secondaryperiod<-function(SPdays=1,rdate=rdata$Date){

  ad1=tapply(rdate,rdata$PPNum,function(x)ceiling(as.numeric(format(x,'%j'))/SPdays))
  ad2=lapply(ad1,function(x)as.numeric(factor(x,ordered=TRUE,levels=unique(x))))
  unlist(ad2)
}


###################
###
### return data (take, effort, effort per area, basis vectors) for all data
###
###################

get_site_data <- function(data){
  require(reshape2)
  require(splines)

  Areaper_r <- Xd_r <- ymat_r <- gbe_r <- tibble()

  properties <- unique(data$Property)

  for(i in seq_along(properties)){

    rdata <- data |>
      filter(Property == properties[i])

    rdata$Date=as.POSIXct(rdata$Date,format="%m/%d/%Y %H:%M")

    ### Create a column for the primary period of interest starting at the beginning of the study, must be a discrete number of months, standard is one month
    # rdata$PPNum=primaryperiodsequence(PPmonths = 3,rdata$Date)

    ### Make methods as factors and order these based on which would occur earlier in the secondary period (usually day)
    rdata$Methods=factor(rdata$Methods,levels=c("Trap","Aerial","Shooting"),ordered = TRUE)

    ### Order the data by site then by date, then by method
    rdata=rdata[order(rdata$Site,rdata$Date,rdata$Methods),]

    ### Set the secondary period (SPdays = number of days for secondary period)
    # rdata$SPNum=secondaryperiod(SPdays=0.1,rdate=rdata$Date)
    rdata <- rdata |>
      group_by(PPNum) |>
      mutate(SPNum = 1:n()) |>
      ungroup()

    ### Upload or create an area data frame giving the area for each site in km2, allow for a buffer if removal efforts on the edge of the study area
    property_area <- rdata$area_km2[1]
    ardf=data.frame(Site=1,Area=property_area)

    ### Set up the area of impact for each method, and method names
    detdf=data.frame(Method=c("Trap","Aerial","Shooting","0"),Area=c(15,40,5,0))
    capnames=detdf$Method[-dim(detdf)[1]]


    ###
    ### Processing the data for input in the MCMC model
    ###
    mptake=dcast(data = rdata[,c("Site","PPNum","Take","Methods","SPNum")],Site+PPNum~SPNum+Methods,value.var = "Take",fun.aggregate = sum)
    mpeff=dcast(data = rdata[,c("Site","PPNum","Effort","Methods","SPNum")],Site+PPNum~SPNum+Methods,value.var = "Effort",fun.aggregate = sum)

    passes=max(apply(mptake,1,function(x)length(which(x>0))))-2
    mths=length(unique(mptake$MYNum))
    sites=max(rdata$Site)

    ### Really slow way of processing the data but it works
    takeinfo=matrix(0,dim(mptake)[1],passes)
    effortinfo=matrix(0,dim(mptake)[1],passes)
    metheff=matrix("0",dim(mptake)[1],passes)

    for(t in 1:dim(mptake)[1]){
      tkind=which(mptake[t,]>0)[-(1:2)]
      takeinfo[t,seq(1,length(tkind),1)]=as.numeric(mptake[t,tkind])
      effortinfo[t,seq(1,length(tkind),1)]=as.numeric(mpeff[t,tkind])
      metheff[t,seq(1,length(tkind),1)]=gsub(".*_","",names(mptake)[tkind])
    }
    ymat=cbind(mptake[,1:2],takeinfo)
    gbe=(effortinfo)
    gtypes=metheff
    areatot=ardf$Area
    sml=dim(ymat)[1]

    ### Use the function to calculate the area of impact, this will work with and without lat/long coordinates
    Areaper=area.of.impact.finder2(as.data.frame(rdata),gtypes,ymat,ardf,gbe,areatot,GPS=TRUE,detdf,plotareas=FALSE)

    ### Create a design matrix for the growth rate analysis, Here I used an intercept (necessary) and basis function expansion on time to allow for growth rate to vary with time
    Xdat=cbind(1,(ymat$PPNum)/max(ymat$PPNum))
    Xd=cbind(1,bs(Xdat[,2],5))

    prop_name <- rdata$Property[1]

    make_long <- function(df){
      df |>
        as_tibble() |>
        mutate(ts = 1:n(),
               Property = prop_name) |>
        pivot_longer(cols = -c(ts, Property))
    }

    Areaper_long <- make_long(Areaper)
    Areaper_r <- bind_rows(Areaper_r, Areaper_long)

    Xd_long <- make_long(Xd)
    Xd_r <- bind_rows(Xd_r, Xd_long)

    ymat_long <- make_long(ymat)
    ymat_r <- bind_rows(ymat_r, ymat_long)

    gbe_long <- make_long(gbe)
    gbe_r <- bind_rows(gbe_r, gbe_long)

  }

  make_wide <- function(df){
    df |>
      pivot_wider(names_from = name,
                  values_from = value) |>
      arrange(Property, ts)
  }

  return(list(
    Areaper = make_wide(Areaper_r),
    Xd = make_wide(Xd_r),
    ymat = make_wide(ymat_r),
    gbe = make_wide(gbe_r)
  ))
}
