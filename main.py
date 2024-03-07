#Commento con il cancelletto
#No commenti multilinea
messaggio = input("Inserisci qualcosa: ")
print(messaggio)

if 1<5:
    print("Ciao")
    #indentazione automatica, mi dice se sono dentro il ciclo o no

#Variabile
#Subito valore da assegnare no scatola vuota a differenza di C++

#x=5
#y=6

#Nomenclatura: regole per definire le variabili
#Illegali sono spazi, trattini, o numeri all'inizio
#Case: pesoPersona, PesoPersona, peso_persona (ultimo suggerito)

citta=["roma","milano","napoli"]
x,y,z= citta

print(x)
print(y)
print(z)

#tipo di dato
x = 5.5
print(type(x)) #stampa il tipo di variabile la funzione type()

x="ciao"
print(type(x)) #stampa il tipo di variabile la funzione type()
#una variabile è di un tipo rispetto al valore che assume (più versatile ma non è semplice riconoscere la tipologia di variabile)

x=False #Lettere maiuscola la prima

#Liste, non sono array ma sono collezione di dati (Array in numpy)

#list: x=["roma","milano","napoli"]
#tuple: x=("roma","milano","napoli")
#range: x=range(6) (a list of number)
#dict: x={"nome":"Luca", "eta": 25}
#set: x={"roma","milano","napoli"}

#si distinguono in base alle proprietà
######################################
#Castin: converzione di un tipo in un altro

x = "5"
y = str(5) #casting

#print(x+y) Errore (concatenazione): il computer vede una stringa
#print(y+x) Errore (sommare intero a una stringa) 

x = 5
y= int("5")
x= float(5)

#######################################
#Stringhe
x="ciao"
y='ciao' #' o " uguale basta non mischiarli

#Per creare stringhe multilinea:
x="""ciao
vado a capo"""
y='''ciao
a capo facile anche con singolo apice'''

#NB stringa è un array di caratteri
print(x[0]) #output: c (anche lo spazio è un carattere)
print(len(x))#len() da la lunghezza della stringa

for carattere in "computer":
    print(carattere)

print(x[:3])
print(x[1:3]) #1 partenza compreso, 7 arrivo escluso
print(x[-4]) #indici negativi 

#Metodi per str
x.upper() #upper case (maiuscolo)
x.lower() #lower case (minuscolo)
x.strip() #toglie lo spazio davanti e in fondo alla stringa
x.replace("o","w") #al posto della 0 ci mette la w

#Concatenare le stringhe
print(x+y) #dove x e y sono due stringhe

#Combinare stringhe e numeri
x = 23
y = 70
z = 1.70
prova = "Ciao sono Luca e sono nato il {2}, peso {0} e altezza {1}" #indici z,x,y
print(prova.forma(x,y,z))

#Escape dei caratteri
prova = "Sono \"figo""
prova = 'Sono alla ricerca dell\'amore'

###############################################

#Boolean

###############################################

#if
x==10
if x<10:
    print("x è minore di 10")
elif x==10:
    print("uguale a 10")
else:
    print("condizione NON verificata")
#and or

###############################################

#while
i=0
while i<6:
    print(i)
    i+=1

i=0
while i<6:
    print(i)
    if i==3:
        break
    i+=1

i=0
while i<6:
    print(i)
    if i==3:
        continue #salta il 3
    i+=1

i=0
while i<6:
    print(i)
    i+=1
else:
    print("ho finito")

###############################################
#for

lista_citta=["milano","roma","napoli"]

for citta in lista_citta:
    print(citta) #come i vector in C++

for x in range(6):
    print(x)

##############################################

#










