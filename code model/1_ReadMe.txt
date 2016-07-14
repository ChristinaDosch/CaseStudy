Zum Testen der verschiedenen Ansätze die start-Funktionen aufrufen:

* start_SO_closed für den SO Ansatz mit geschlossener Formel für den Erwartungswert
* start_SO_discretization für den SO Ansatz mit Diskretisierung, OHNE Batterie
* start_RO für den RO Ansatz OHNE Batterie
* start_ROvsSO für einen Vergleich von RO ohne Batterie und SO closed-form ohne Batterie
* start_SO_discr_bat für den SO Ansatz mit brute force approach und Diskretisierung (brute force)
* start_SO_discr_bat_smart_constr für den SO Ansatz mit smart approach und Diskretisierung ( O(K*T) Variablen )
* start_SO_discr_bat_smart_constr für den SO Ansatz mit smart approach und Diskretisierung ( inneres Optimierungsproblem )

In init_constraints werden die Constraints definiert, 
in init_parameters legen wir die Parameter fest. 