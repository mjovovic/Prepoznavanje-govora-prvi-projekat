Prepoznavanje govora prvi projekat

Alat koji prikazuje frekventni spektar unapred snimljene reči.


Ulaz
Alat mora na ulazu da obezbedi sledeće stavke:
Korišćenje sirovog .wav fajla kao ulazni zvučni signal. Na snimku se nalazi izgovorena jedna reč, s time da je moguće da pre i posle reči postoji tišina proizvoljne dužine. Ovu tišinu treba iseći pre početka obrade, primenom algoritma opisanog na vežbi 2. Sistem treba da prijavi poruku o grešci ako smatra da na snimku nema izgovorene reči.
Izbor Trouglaste, Hamming, Hanning, Blackman, ili ni jedne prozorske funkcije.
Izbor širine DFT prozora.
Obezbediti makar pet .wav datoteka za testiranje, svaka dužine po par sekundi:
Govorni signal - muški glas.
Govorni signal - ženski glas.
Jedan ton na muzičkom instrumentu.
Nasumični šum (za proveru poruke o grešci).
Generisan signal koji sadrži nekoliko probranih harmonika (napisati jednostavan program za generisanje ovog signala).
Izlaz
Alat mora da omogući prikaz frekventnog spektra na dva načina:
Prikaz čitavog signala iz .wav fajla i označena lokacija početka i kraja reči.
Prikaz frekventnog spektra za jedan prozor, gde se spektar prikazuje kao histogram.
Prikaz frekventnog spektra čitavog signala, gde se spektar prikazuje kao sonogram.
