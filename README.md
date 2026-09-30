# Simulador de propulsión D–³He

Simulador educativo que compara dos naves con datos publicados:

- **Nave A · clase Dawn**: propulsión iónica de xenón. Voló entre 2007 y 2018.
- **Nave B · Direct Fusion Drive de 2 MW**: propulsión por fusión deuterio–helio-3. Existe solo como diseño.

Proyecto de la Feria Científica Escolar 2026, Kingston College (Concepción, Chile). Acompaña al paper «Propulsión por fusión D–³He frente a propulsión iónica de xenón».

**Probarlo:** abre `index.html` en cualquier navegador, o la versión publicada en GitHub Pages.

## Cómo está construido

El simulador tiene dos capas separadas:

1. **Capa física (determinista).** El objeto `Fisica`, al inicio del script de `index.html`, ejecuta solo tres ecuaciones cerradas:
   - Ecuación del cohete: Δv = vₑ · ln(m₀ / m_f)
   - Empuje y potencia: F = 2ηP / vₑ
   - Tiempo de encendido: t_b = (m₀ vₑ / F) · [1 − e^(−Δv/vₑ)]
2. **Capa de interfaz (juego).** Misiones, cinemáticas, estrellas y gráficos. No calcula física: solo muestra lo que entrega la capa física. Las cinemáticas son ilustrativas y no están a escala.

No usa inteligencia artificial en tiempo de ejecución. La pista de la animación mide Δv ganado, no distancia.

## Validación

La pestaña **Validación** ejecuta 10 casos de prueba y los compara con el modelo en Python (`modelo/modelo_propulsion_he3.py`, resultados en `modelo/results.json`) y con los valores publicados. Criterio: error relativo menor que 1 %.

| Caso | Publicado | Simulador | Error |
|---|---|---|---|
| DFD a Haumea, Δv | 64,16 km/s (Aime et al., 2021) | 64,15 km/s | 0,02 % |
| DFD a Haumea, tiempo de empuje | 510 días | 510,1 días | 0,01 % |
| Motor NEXT, eficiencia | 0,71 (Herman et al., 2008) | 0,706 | 0,61 % |

## Fuentes principales (acceso gratuito)

- Aime, Gajeri y Kezerashvili (2021). *Exploration of trans-Neptunian objects using the Direct Fusion Drive.* Acta Astronautica 178. https://arxiv.org/abs/2009.12633
- Rayman et al. (2006). *Dawn: A mission in development.* https://science.nasa.gov/wp-content/uploads/2023/09/Dawn_overview.pdf
- Herman, Soulas y Patterson (2008). *NEXT ion thruster.* https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20080012739.pdf
- Galea et al. (2023). *The Princeton Field-Reversed Configuration.* https://w3.pppl.gov/ppst/docs/galea2023jfe.pdf

## Límites

- No calcula trayectorias, gravedad ni destinos.
- Supone empuje e impulso específico constantes.
- No modela el reactor de fusión: la dificultad de ignición viene del cálculo de Lawson del paper.
