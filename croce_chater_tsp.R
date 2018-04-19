
#install.packages("geosphere")
#install.packages("mapdata")
#install.packages("maps")


library(gtools)

######################################################################
#Algorithme montrant le  temps de recherche selon le nombre de ville #
######################################################################

graphic_temps_recherche_selon_nb_ville<-function(nb_ville){
  #Data.frame de stockage des temps de calcul selon le nombre de villes
  prog_temps<-data.frame(matrix(NA, nrow = nb_ville, ncol = 2))
  colnames(prog_temps)<-c("Nb ville","Tps computing")
  for(v in 1:nb_ville){
    #Modelisation des villes par 2 lois unif
    #Nombre de ville
    m=v
    #Creation des villes sur une carre 1*1 par 2 loi unif
    x1=runif(m, min=0,max=1)
    x2=runif(m, min=0,max=1)
    ville=cbind(x1,x2)
    
    #Calcul de la matrice des distances entre les villes
    D=dist(ville)
    D=as.matrix(D)
    
    #Démarage du chronomètre
    start.time <- Sys.time()
    #debut de la fonction de recherche exhaustive
    recherche_plus_court_chemin(m,D)
    #fin du chronomètre
    end.time <- Sys.time()
    prog_temps[v,1]<-v
    prog_temps[v,2]<-as.numeric(as.numeric(end.time)-as.numeric(start.time))
    
  }
  plot(x=prog_temps[,1],y=prog_temps[,2])
  lines(x=prog_temps[,1],y=prog_temps[,2], type = "l",col="red", 
        main="Temps de calcul selon le nombre de villes",xlab="Nombre de ville",ylab="Durée en secondes")
}

##############################################
# Fonction de calcul de distance d'un chemin #
##############################################
## chemin_test=chemin a calculer, m=nb ville, D=Matrice des distances
distance<-function(chemin_test,m,D){
  tps<-0
  for(i in 1:m){
    tps<-tps+D[chemin_test[i],chemin_test[i+1]]
  }
  return<-tps
}


############################
# Algorithme recuit_simulé #
############################

#Algorithme recuit simulé
library(car)
library(gtools)
library(sp)
library(geosphere)





#Algorithme de Metropolis
metropolis<-function(d1, c1, d2, c2, T){
  if(d1 <= d2){
    c2 <- c1
    d2 <- d1
  } else {
    if(runif(1,0,1) > exp(-(d1-d2)/T)){
      c2 <- c1
      d2 <- d1
    }
  }
  return(c2)
}
  

#########################
# Algo du recuit simule #
#########################

## (m= nb ville, D= matrice des distances, nbiteration = nombre d'iteration de l'algorithme)
recuit_simule<-function(m,D, nbiteration){
  #x0
  chemin_act<-c(seq(1,m,1),1)
  boucle<-log(nbiteration)/100
  d_act<-distance(chemin_act,m,D)
  #Inversement possible
  val_inverser<<-combinations((m),2,chemin_act)
  #iteration des k
  for(k in 1:100){
    #Temperature de cette sequence (varie)
    T<-1/k
    for(n in ceiling(exp((k-1)*boucle)):(ceiling((exp(k*boucle)))+1)){ #T est constante par morceaux
      #numero du chemin a tester
      num_test<<-sample(1:nrow(val_inverser),1)
      
      #Selection du chemin voisin a tester
      chemin_test<-recode(chemin_act, 'val_inverser[num_test,1]=val_inverser[num_test,2]; val_inverser[num_test,2]=val_inverser[num_test,1]')
      
      #Distance du chemin actuel
      d_test<-distance(chemin_test,m,D)
      
      #Probabilite de sauter
      chemin_act <- metropolis(d_act,chemin_act, d_test, chemin_test)
      d_act<-distance(chemin_act,m,D)
    }
  }
  return<-list(chemin_court=chemin_act, distance=distance(chemin_act,m,D))
}

############################################################################
# Représentation graphique du parcours des états à l'aide du recuit simulé #
############################################################################

# m= nb ville,D= matrice des distances, nbiteration=nombre d'iteration de l'algorithme
recuit_simule_graphique<-function(m,D, nbiteration){
  
  #x0
  chemin_act<-c(seq(1,m,1),1)
  boucle<-log(nbiteration)/100
  d1<-distance(chemin_act,m,D)
  resultat<-chemin_act
  #creation plot
  plot(x=0,y=d1,xlim=c(1,nbiteration), ylim=c(d1/2,(d1*1.25)), xlab="Numéro de l'itération", ylab="Distance du chemin à la fin de l'itération n")
  
  #Invertissement possible
  val_inverser<<-combinations((m),2,chemin_act)
  indice_graph<-1
  #iteration des k
  for(k in 1:100){
    #Temperature de cette sequence
    T<-1/k
    for(n in ceiling(exp((k-1)*boucle)):(ceiling((exp(k*boucle)))+1)){
      #numero du chemin a tester
      num_test<<-sample(1:nrow(val_inverser),1)
      #Selection du chemin voisin a tester
      chemin_test<-recode(chemin_act, 'val_inverser[num_test,1]=val_inverser[num_test,2]; val_inverser[num_test,2]=val_inverser[num_test,1]')
      
      #Distance du chemin actuel
      d2<-distance(chemin_test,m,D)
      #Probabilite de sauter
      p<-exp((1/T)*(d1-d2))
      if(p>=1){
        chemin_act<-chemin_test
        d1<-d2
      }else{
        if(runif(1,0,1)<p){
          chemin_act<-chemin_test
          d1<-d2
        }
      }
      if(distance(resultat,m,D)>d2){
        resultat<-chemin_test
      }
      points(x=indice_graph,y=d1)
      indice_graph<-indice_graph+1
      
    }
    
  }
  return<-list(chemin_court=resultat, distance=distance(resultat,m,D))
}



##################################
# Recuit_simulé_temperature_fixe #
##################################


#####################################
# Algo du recuit a temperature fixe #
#####################################
# m= nb ville,D= matrice des distances, nbiteration=nombre d'iteration de l'algorithme, temperature=température
temperature_fixe<-function(m,D, nbiteration,temperature){
  #x0
  chemin_act<-c(3, 11,  2, 19, 12,  6,  5,  7,  8, 10,  9,  1, 20, 21, 22, 16, 17, 18,  4, 15, 13, 14,  3)
  boucle<-log(nbiteration)/100
  d1<-distance(chemin_act,m,D)
  resultat<-chemin_act
  
  #Invertissement possible
  val_inverser<<-combinations((m),2,chemin_act)
  indice_graph<-1
  #iteration des k
  for(k in 1:100){
    #Temperature de cette sequence
    
    for(n in ceiling(exp((k-1)*boucle)):(ceiling((exp(k*boucle)))+1)){
      #numero du chemin a tester
      num_test<<-sample(1:nrow(val_inverser),1)
      #Selection du chemin voisin a tester
      chemin_test<-recode(chemin_act, 'val_inverser[num_test,1]=val_inverser[num_test,2]; val_inverser[num_test,2]=val_inverser[num_test,1]')
      
      #Distance du chemin actuel
      d2<-distance(chemin_test,m,D)
      #Probabilite de sauter
      p<-exp((1/temperature)*(d1-d2))
      if(p>=1){
        chemin_act<-chemin_test
        d1<-d2
      }else{
        if(runif(1,0,1)<p){
          chemin_act<-chemin_test
          d1<-d2
        }
      }
      if(distance(resultat,m,D)>d2){
        resultat<-chemin_test
      }
      
    }
    
  }
  return<-list(chemin_court=resultat, distance=distance(resultat,m,D))
}
  
#################################################################################################
# Représentation graphique du parcours des états avec l'algorithme du recuit a température fixe #
#################################################################################################
# m= nb ville,D= matrice des distances, nbiteration=nombre d'iteration de l'algorithme, temperature=temperature fixe
temperature_fixe_graphique<-function(m,D, nbiteration,temperature){
  
  #x0
  chemin_act<-c(seq(1,m,1),1)
  boucle<-log(nbiteration)/100
  d1<-distance(chemin_act,m,D)
  resultat<-chemin_act
  #creation plot
  plot(x=0,y=d1,xlim=c(1,nbiteration), ylim=c(d1/2,(d1*1.25)), xlab="Numéro de l'itération", ylab="Distance du chemin à la fin de l'itération n")
  
  #Invertissement possible
  val_inverser<<-combinations((m),2,chemin_act)
  indice_graph<-1
  #iteration des k
  for(k in 1:100){
    #Temperature de cette sequence
    T<-1/k
    for(n in ceiling(exp((k-1)*boucle)):(ceiling((exp(k*boucle)))+1)){
      #numero du chemin a tester
      num_test<<-sample(1:nrow(val_inverser),1)
      #Selection du chemin voisin a tester
      chemin_test<-recode(chemin_act, 'val_inverser[num_test,1]=val_inverser[num_test,2]; val_inverser[num_test,2]=val_inverser[num_test,1]')
      
      #Distance du chemin actuel
      d2<-distance(chemin_test,m,D)
      #Probabilite de sauter
      p<-exp((1/temperature)*(d1-d2))
      if(p>=1){
        chemin_act<-chemin_test
        d1<-d2
      }else{
        if(runif(1,0,1)<p){
          chemin_act<-chemin_test
          d1<-d2
        }
      }
      if(distance(resultat,m,D)>d2){
        resultat<-chemin_test
      }
      points(x=indice_graph,y=d1)
      indice_graph<-indice_graph+1
      
    }
    
  }
  return<-list(chemin_court=resultat, distance=distance(resultat,m,D))
}

#################################################################################################

########################################################
# Question 1 : Recherche exhaustive du meilleur chemin #
########################################################

#Modélisation des villes par 2 lois unif
#Nombre de ville
m=5
#Création des villes sur une carré 1*1 par 2 loi unif
x1=runif(m, min=0,max=1)
x2=runif(m, min=0,max=1)
ville=cbind(x1,x2)

#Calcul de la matrice des distances entre les villes
D=dist(ville)
D=as.matrix(D)




#### Fonction de recherche exhaustive du plus court chemin ####
## retourne le chemin le plus court et sa distance   
## (m=nb_ville,D=matrice des distance)
recherche_plus_court_chemin<-function(m,D){
  #Recherche de tout les chemins possible
  x=seq(1,m,1)
  chemin=permutations(m,m,x) ##Diviser par deux ?
  chemin=cbind(chemin,chemin[,1])
  #Calcul du temps des chemins
  distance<-rep(NA,nrow(chemin))
  for(j in 1:nrow(chemin)){
    tps<-0
    for(i in 1:m){
      tps<-tps+D[chemin[j,i],chemin[j,i+1]]
    }
    distance[j]<-tps
  }
  #déteminer le chemin le plus court
  return<-list(chemin_court=chemin[which.min(distance),],distance=min(distance))
}


### Algorithme montrant le  temps de recherche selon le nombre de ville ###
graphic_temps_recherche_selon_nb_ville<-function(nb_ville){
  #Data.frame de stockage des temps de calcul selon le nombre de villes
  prog_temps<-data.frame(matrix(NA, nrow = nb_ville, ncol = 2))
  colnames(prog_temps)<-c("Nb ville","Tps computing")
  for(v in 1:nb_ville){
    #Modelisation des villes par 2 lois unif
    #Nombre de ville
    m=v
    #Creation des villes sur une carre 1*1 par 2 loi unif
    x1=runif(m, min=0,max=1)
    x2=runif(m, min=0,max=1)
    ville=cbind(x1,x2)
    
    #Calcul de la matrice des distances entre les villes
    D=dist(ville)
    D=as.matrix(D)
    
    #Démarage du chronomètre
    start.time <- Sys.time()
    #debut de la fonction de recherche exhaustive
    recherche_plus_court_chemin(m,D)
    #fin du chronomètre
    end.time <- Sys.time()
    prog_temps[v,1]<-v
    prog_temps[v,2]<-as.numeric(as.numeric(end.time)-as.numeric(start.time))
    
  }
  plot(x=prog_temps[,1],y=prog_temps[,2])
  lines(x=prog_temps[,1],y=prog_temps[,2], type = "l",col="red", 
        main="Temps de calcul selon le nombre de villes",xlab="Nombre de ville",ylab="Durée en secondes")
}


#Lancer la fonction  
a<-recherche_plus_court_chemin(m,D)
#resultat
a$chemin_court
a$distance

#graphique du temps de recherche selon le nombre de ville
graphic_temps_recherche_selon_nb_ville(m)



#Nombre de ville
#Position des villes
data<-read.table(file="/home/gis3/cchater/Linux/GIS4/PS/20villes.txt", header=T)

##############
# Question 2 #
##############

#http://alain.camanes.free.fr/universite/vulgarisation/fs07/voyageur_complement.pdf

#########################################################################
# Question 3 : Algorithme du recuit simulé pour le voyageur de commerce #
#########################################################################


#Comparaison de la recherche du meilleur chemin avec la méthode exhaustive et le recuit simulé pour m petit

#Modelisation des villes par 2 lois unif
#Nombre de ville
m=3
#Creation des villes sur une carre 1*1 par 2 loi unif
x1=runif(m, min=0,max=1)
x2=runif(m, min=0,max=1)
ville=cbind(x1,x2)

#Calcul de la matrice des distances entre les villes
D=dist(ville)
D=as.matrix(D)
plot(ville)
D

#Lancer la recherche exhaustive
a<-recherche_plus_court_chemin(m,D)
#resultat
a$chemin_court
a$distance

#Lancer la recherche du plus court chemin par la methode mcmc
b<-recuit_simule(m,D, 5000)
#resultat
b$chemin_court
b$distance

#Recherche du meilleur chemin avec l'algorithme du recuit_simulé
#Matrice des distance
m<-22
D=distm(data[,3:2])
D=as.matrix(D)

chemin_court<-recuit_simule(m, D, 10000)
chemin_court$chemin_court
chemin_court$distance

#Représentation graphique du parcours des états à l'aide du recuit simulé
#Matrice des distance

chemin_court<-recuit_simule_graphique(m, D, 5000)



#Recherche du meilleur chemin avec l'algorithme du recuit a température fixe 
#Matrice des distance
D=distm(data[,3:2])
D=as.matrix(D)
D<-D/max(D)


chemin_court<-temperature_fixe(m, D, 100000, 0.05)


#Représentation graphique du parcours des états avec l'algorithme du recuit a température fixe 
#Matrice des distance
D=distm(data[,3:2])
D=as.matrix(D)
D<-D/max(D)
chemin_court<-temperature_fixe_graphique(m, D, 5000, 0.05)


##############
# Question 4 #
##############
#Nombre de ville
library(maps)
library(mapdata)
m<-22
#Position des villes
data<-read.table(file="/home/gis3/cchater/Linux/GIS4/PS/20villes.txt", header=T)
D=distm(data[,3:2])
D=as.matrix(D)
D_max<-max(D)
D<-D/D_max
chemin_court<-temperature_fixe(m,D, 1000,0.05)
D<-D*D_max
chemin_court$distance<-distance(chemin_court$chemin_court,m,D)
chemin_court$distance
a<-chemin_court
#Modélisation de la map

map('france', border=0, fill=TRUE, 
    bg="gray", col="white",mar=rep(0,4))
points(data[,3],data[,2],pch=20,col="blue")
text(data[,3], data[,2], data[,1], cex=0.6, pos=1, col="red")
map.scale()

for(i in 1:(m)){
  
  segments(data[a$chemin_court[i],3], data[a$chemin_court[i],2], data[a$chemin_court[i+1],3],data[a$chemin_court[i+1],2])
  
}
