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

#castin: converzione di un tipo in un altro

x = "5"
y = str(5) #casting

#print(x+y) Errore (concatenazione): il computer vede una stringa
#print(y+x) Errore (sommare intero a una stringa) 

x = 5
y= int("5")
x= float(5)