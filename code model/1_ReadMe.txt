Zum Testen der verschiedenen Ans�tze die start-Funktionen aufrufen:

* start_SO_closed f�r den SO Ansatz mit geschlossener Formel f�r den Erwartungswert
* start_SO_discretization f�r den SO Ansatz mit Diskretisierung, OHNE Batterie
* start_RO f�r den RO Ansatz OHNE Batterie
* start_ROvsSO f�r einen Vergleich von RO ohne Batterie und SO closed-form ohne Batterie
* start_SO_discr_bat f�r den SO Ansatz mit brute force approach und Diskretisierung (brute force)
* start_SO_discr_bat_smart_constr f�r den SO Ansatz mit smart approach und Diskretisierung ( O(K*T) Variablen )
* start_SO_discr_bat_smart_constr f�r den SO Ansatz mit smart approach und Diskretisierung ( inneres Optimierungsproblem )

In init_constraints werden die Constraints definiert, 
in init_parameters legen wir die Parameter fest. 