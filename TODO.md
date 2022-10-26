# TODO

## Hinweise für die Abgabe
- [ ] in jedem Zeitschritt vti schreiben
- [ ] auf Terminal in jedem Zeitschritt infos ausgeben, z.B. residual
- [ ] 1 o. 2 Szenarien gibt es zum testen, 1 o. 2 Szenarien unbekannt - [ ]> nur Ergebnis
- [ ] normal residual 9e- [ ]6, 60- [ ]25 Solver Schritte, 383 time steps
- [ ] Solver Schritte nehmen kontinuierlich ab (bis auf letzten)
- [ ] last line: time step 383, t: 9.98602/10, dt: 0.0139948, res. 
- [ ] vti files in jedem Zeitschritt im out Ordner
- [ ] in paraview anschauen
- [ ] wenn alles passt: http://opendihu- [ ]ci.informatik.uni- [ ]stuttgart.de/numsim/
- [ ] Nach Upload manuell Seite mit F5 refreshen - [ ]> dann output sehen
- [ ] scp - [ ]r ... ipvslogin: ...
- [ ] kommt ins ipvslogin /home, aber ist überall im cluster auf /home sichtbar
- [ ] auf sgscl1 kann kompliert werden, aber besser nichts ausführen
- [ ] slurm nutzen um auf eigentliche knoten des clusters zu kommen
- [ ] (script.sh erstellen mit "srun hostname\n sleep 1000000")
- [ ] dann sbatch script.sh
- [ ] aktuelle jobs ansehen: squeue, output ist in slurm- [ ]1345.out o.ä.
- [ ] job anhalten: scancel 1345
- [ ] von website script script_parallel.sh holen
- [ ] interaktiv ausführen: srun - [ ]- [ ]pty bash
