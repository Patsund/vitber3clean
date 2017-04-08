import oppgave1 as o1
import oppgave2 as o2
import oppgave3 as o3

def printOptions():
	outStr = "\n1) Oppgave 1"
	outStr+= "\n2) Oppgave 2"
	outStr+= "\n3) Oppgave 3"
	outStr+= "\n4) Velg mellom å vise og lagre figurer"
	outStr+= "\n0) Avslutt"
	outStr+= "\nDitt valg:"
	return outStr

def printSave():
	outStr = "\n1) Vis figurer"
	outStr+= "\n2) Lagre figurer"
	outStr+= "\n0) Avslutt"
	outStr+= "\nDitt valg:"
	return outStr

if __name__ == "__main__":
	switch = -1
	save = -1
	while(save):
		try:
			save = int(input(printSave()))
		except:
			pass
		if save != 1 and save != 2 and save != 0:
			print("Du oppga ikke en gyldig verdi for å lagre/vise figurer")
		if save == 0:
			switch=False
			break
		else:
			break
	while(switch):
		try:
			switch = int(input(printOptions()))
		except:
			pass
		if save == 3:
			pass
		if switch == 1:
			if save == 1:
				o1.oppgave1()
			elif save == 2:
				o1.oppgave1(True)
		elif switch == 2:
			datapath = 'C:/Users/Patrik/Downloads/NorKyst-800m.nc'
			if save == 1:
				o2.oppgave2(datapath)
			elif save == 2:
				o2.oppgave2(datapath, True)
		elif switch == 3:
			datapath = 'C:/Users/Patrik/Downloads/NorKyst-800m.nc'
			if save == 1:
				o3.oppgave3(datapath)
			elif save == 2:
				o3.oppgave3(datapath, True)
		elif switch == 4:
			while(save):
				try:
					save = int(input(printSave()))
				except:
					pass
				if save != 1 and save != 2 and save != 0:
					print("Du oppga ikke en gyldig verdi for å lagre/vise figurer")
				if save == 0:
					switch=False
					break
				else:
					break
		elif switch == 0:
			break
		else:
			print("Du oppga ikke en gyldig verdi for å velge oppgave\n")