import numpy as np 

# fonction ou le code a reutiliser se trouve
def hello_world (data):
    print ("La valeur recu par notre fonction hello world -   {}".format(data))
    return 

# donnees que je dois generer pour cet exemple
data = np.linspace(0,100)
m = len(data)

# appel de notre fonction hello_world qui effectu le meme travaille a plusieurs reprises
for i in range(m):
    hello_world(data[i])
print ('Fin du message')