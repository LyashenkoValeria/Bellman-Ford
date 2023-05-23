package org.example;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

public class LinearBF {
    public static int RUNS = 10;
    public static int NODE_SIZE = 1500;

    public static class Edge {
        int u;
        int[] v;
        int[] w;
        public Edge(int u, int[] v, int[] w){
            this.u = u;
            this.v = v;
            this.w = w;
        }
    }

    public static Edge[] edges;

    public static void main(String[] args) {
        long sumTime = 0;
        long[] time = new long[RUNS];

        for (int r = 1; r <= RUNS; r++) {
            int[][] graph = graphGeneration();
            getData(graph);

            long startTime1 = System.currentTimeMillis();
            bellmanFord(0);
            long endTime1 = System.currentTimeMillis();

            time[r - 1] = endTime1 - startTime1;
            sumTime += time[r - 1];
            System.out.println("Run - " + r + " | Time: " + time[r - 1]);
        }

        //Расчёт значений
        double mean = ((double) sumTime / (double) RUNS);
        double dispersion = 0;
        System.out.println("Mean time: " + mean);
        for (int r = 0; r < RUNS; r++) {
            dispersion += pow(time[r] - mean, 2);
        }
        dispersion = dispersion/(double) RUNS;
        System.out.println("Dispersion time: " + dispersion);
        double s = sqrt((RUNS*dispersion)/(RUNS-1));
        double radius = 2.26*s/sqrt(RUNS);
        System.out.println("Interval time: [" + (mean-radius) + "; " + (mean+radius) + "]");
    }

    public static void bellmanFord(int src) {
        int[] weights = new int[NODE_SIZE];
        Arrays.fill(weights, (int) 1e7);
        weights[src] = 0;

        for (int i = 0; i < NODE_SIZE; i++) {
            for (int u = 0; u < edges.length; u++) {
                for (int v = 0; v < edges[u].v.length; v++) {
                    if (weights[edges[u].u] + edges[u].w[v] < weights[edges[u].v[v]]) {
                        weights[edges[u].v[v]] = weights[edges[u].u] + edges[u].w[v];
                    }
                }
            }
        }

        for (int i = 0; i < NODE_SIZE; i++) {
            for (int u = 0; u < edges.length; u++) {
                for (int v = 0; v < edges[u].v.length; v++) {
                    if (weights[edges[u].u] + edges[u].w[v] < weights[edges[u].v[v]]) {
                        // System.out.println("Graph contains negative cycle");
                        return;
                    }
                }
            }
        }
    }

    public static int[][] graphGeneration(){
        Random rand = new Random();
        int[][] graph = new int[NODE_SIZE][NODE_SIZE];

        for (int i = 0; i < NODE_SIZE; i++){
            for (int j = 0; j < NODE_SIZE; j++){
                if (i!=j) graph[i][j] = rand.nextInt(0,10001);
            }
        }
        return graph;
    }

    public static void getData(int[][] graph) {
        List<Edge> listOfEdges = new ArrayList<>();
        for (int i = 0; i < NODE_SIZE; i++) {

            List<Integer> listOfV = new ArrayList<>();
            List<Integer> listOfW = new ArrayList<>();
            for (int j = 0; j < NODE_SIZE; j++) {
                if (graph[i][j] != 0) {
                    listOfV.add(j);
                    listOfW.add(graph[i][j]);
                }
            }
            Edge edge = new Edge(i, listOfV.stream().mapToInt(Integer::intValue).toArray(), listOfW.stream().mapToInt(Integer::intValue).toArray());
            listOfEdges.add(edge);
        }

        edges = listOfEdges.toArray(new Edge[NODE_SIZE]);
    }
}