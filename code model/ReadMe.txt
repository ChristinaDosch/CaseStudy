Zum Testen der verschiedenen Ansätze die start-Funktionen aufrufen:

* start_SO_closed für den SO Ansatz mit geschlossener Formel für den Erwartungswert
* start_SO_discretization für den SO Ansatz mit Diskretisierung des Erwartungswert, OHNE Batterie
* start_SO_discretization_battery für den SO Ansatz mit Diskretisierung des Erwartungswert MIT Batterie (brute force)
* start_RO für den RO Ansatz OHNE Batterie
* start_ROvsSO für einen Vergleich von RO ohne Batterie und SO closed-form ohne Batterie

In init_constraints werden die Constraints definiert, 
in init_parameters legen wir die Parameter fest. Hieran kann und sollte man schrauben, je nachdem, welchen Ansatz man verwendet.
Insbesondere die Anzahl T an time steps, die noch zu einer sinnvollen Laufzeit führt, ist in init_parameters für jeden Ansatz vermerkt. 
Bitte zum Testen entsprechend auswählen.