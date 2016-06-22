Zum Testen der verschiedenen Ans�tze die start-Funktionen aufrufen:

* start_SO_closed f�r den SO Ansatz mit geschlossener Formel f�r den Erwartungswert
* start_SO_discretization f�r den SO Ansatz mit Diskretisierung des Erwartungswert, OHNE Batterie
* start_SO_discretization_battery f�r den SO Ansatz mit Diskretisierung des Erwartungswert MIT Batterie (brute force)
* start_RO f�r den RO Ansatz OHNE Batterie
* start_ROvsSO f�r einen Vergleich von RO ohne Batterie und SO closed-form ohne Batterie

In init_constraints werden die Constraints definiert, 
in init_parameters legen wir die Parameter fest. Hieran kann und sollte man schrauben, je nachdem, welchen Ansatz man verwendet.
Insbesondere die Anzahl T an time steps, die noch zu einer sinnvollen Laufzeit f�hrt, ist in init_parameters f�r jeden Ansatz vermerkt. 
Bitte zum Testen entsprechend ausw�hlen.