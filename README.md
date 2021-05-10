# Identify COVID-19 Hidden Spreaders over Contact Tracing Networks

This python3 code is a demo for the paper **Identify COVID-19 Hidden Spreaders over Contact Tracing Networks**. In this demo, a spreading of COVID-19 can be simulated on an real-world network "[INFECTIOUS: STAY AWAY](http://konect.cc/networks/sociopatterns-infectious/)". Based on the symptomatic infections, the infection risk of other nodes can be ranked using the proposed method.

## Install required packages

```
pip install -r requirements.txt
```

It takes minutes to install.

## Simulate COVID-19 spreading

```
python generate.py
```

The spreaing parameters can be changed in `model.py`.

## Rank potential infections

```
python rank.py
```

It takes minutes to process. The result can be shown in `out.txt` where the first column is the index of nodes, the second column is the theoretical possibility of infection, the third column is a bool value to show if the node is infected in the simulation. In addition, a visualization of the network shows the risk of infection from low-risk nodes (in green) to high-risk nodes (in red).
