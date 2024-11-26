Top->NA        Third->Fifth   Fifth->NA      Fifth->Third   Third->NA      Second->Second Second->NA
 [8] Fifth->Second  Sixth->Fifth   Fourth->Top    Second->Fifth  Sixth->Sixth   Sixth->NA      Fifth->Top
[15] Top->Sixth     Top->Top       Fourth->Sixth  Second->Top    Third->Top     Top->Fourth    Fourth->NA
[22] Third->Third   Third->Second  Fourth->Fifth  Fifth->Fifth   Fourth->Second Second->Fourth Second->Third
[29] Third->Fourth  Sixth->Fourth  Fourth->Third  Top->Fifth     Fourth->Fourth Fifth->Sixth   Third->Sixth
[36] Second->Sixth  Top->Second    Sixth->Third   Top->Third     Fifth->Fourth  Sixth->Second  Sixth->Top
42 Levels: Top->NA Third->Fifth Fifth->NA Fifth->Third Third->NA Second->Second Second->NA ... Sixth->Top

[1] Low->        High->       Mid->
[11]


"High->High", "High->Mid", "High->Low", "Mid->High", "Mid->Mid", "Mid->Low", "Low->High", "Low->Mid", "Low->Low"


user_input = input('Name: ')

if user_input:
	name = user_input
else:
	name = 'N/A'

print(name)

#is eq to

user_input = input('Name: ')
name = user_input if user_input else 'N/A'

#also eq to
name = user_input or 'N/A'