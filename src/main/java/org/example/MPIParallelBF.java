package org.example;

import mpi.MPI;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import static mpi.MPI.*;

public class MPIParallelBF {
    public static int RUNS = 10;
    public static int NODE_SIZE = 1500;

    static ThreadData threadData;
    static boolean negativeCycle = false;
    public static int rank, size;

    public static class Edge {
        int u;
        int[] v;
        int[] w;

        public Edge(int u, int[] v, int[] w) {
            this.u = u;
            this.v = v;
            this.w = w;
        }
    }

    public static class ThreadData {
        Edge[] edges;
        int[] weights;
        int[] additional;

        public ThreadData(Edge[] edges, int[] weights, int[] additional) {
            this.edges = edges;
            this.weights = weights;
            this.additional = additional;
        }
    }

    public static void main(String[] args) {
        MPI.Init(args);

        long sumTime = 0;
        long[] time = new long[RUNS];

        for (int r = 1; r <= RUNS; r++) {
            rank = MPI.COMM_WORLD.Rank();
            size = MPI.COMM_WORLD.Size();
            int[][] graph = new int[NODE_SIZE][NODE_SIZE];

            if (rank == 0){
                graph = graphGeneration();
            }
            for (int i=0; i < graph.length; i++){
                MPI.COMM_WORLD.Bcast(graph[i], 0, graph[i].length, INT, 0);
            }
            getData(graph);

            long startTime1 = 0;
            if (rank == 0) {
                startTime1 = System.currentTimeMillis();
            }

            bellmanFord(0);

            if (rank == 0) {
                long endTime1 = System.currentTimeMillis();
                time[r - 1] = endTime1 - startTime1;
                sumTime += time[r - 1];
                System.out.println("Run - " + r + " | Time: " + time[r - 1]);
            }
        }

        if (rank == 0) {
            //Расчёт значений
            double mean = ((double) sumTime / (double) RUNS);
            double dispersion = 0;
            System.out.println("Mean time: " + mean);
            for (int r = 0; r < RUNS; r++) {
                dispersion += pow(time[r] - mean, 2);
            }
            dispersion = dispersion / (double) RUNS;
            System.out.println("Dispersion time: " + dispersion);
            double s = sqrt((RUNS * dispersion) / (RUNS - 1));
            double radius = 2.26 * s / sqrt(RUNS);
            System.out.println("Interval time: [" + (mean - radius) + "; " + (mean + radius) + "]");
        }
        MPI.Finalize();
    }

    public static int[][] graphGeneration() {
        Random rand = new Random();
        int[][] graph = new int[NODE_SIZE][NODE_SIZE];

        for (int i = 0; i < NODE_SIZE; i++) {
            for (int j = 0; j < NODE_SIZE; j++) {
                if (i != j) graph[i][j] = rand.nextInt(0, 10001);
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
        threadData = new ThreadData(listOfEdges.toArray(new Edge[NODE_SIZE]), new int[NODE_SIZE], new int[NODE_SIZE]);
    }

    public static void bellmanFord(int src) {
        int step = NODE_SIZE / size;
        int firstNode = rank * step;
        int lastNode = rank == size - 1 ? NODE_SIZE : firstNode + step;

        for (int i = 0; i < NODE_SIZE; i++) Arrays.fill(threadData.weights, (int) 1e7);
        threadData.weights[src] = 0;
        threadData.additional = threadData.weights;

        for (int i = 0; i < NODE_SIZE+1; i++) {
            negativeCycle = i == NODE_SIZE;
            threadBellmanFord(firstNode, lastNode);
            MPI.COMM_WORLD.Reduce(threadData.additional,0, threadData.weights, 0,
                    threadData.additional.length, INT, MIN, 0);
            MPI.COMM_WORLD.Bcast(threadData.weights, 0, threadData.weights.length, INT, 0);
        }
    }


    public static void threadBellmanFord(int firstNode, int lastNode) {
        for (int start = firstNode; start < lastNode; start++) {
            for (int end = 0; end < threadData.edges[start].v.length; end++) {
                if ((threadData.weights[threadData.edges[start].u]
                        + threadData.edges[start].w[end])
                        < threadData.weights[threadData.edges[start].v[end]]) {
                    if (negativeCycle) {
                        // System.out.println("Graph contains negative cycle");
                    } else {

                        int res = threadData.weights[threadData.edges[start].u]
                                + threadData.edges[start].w[end];
                        threadData.additional[threadData.edges[start].v[end]] = res;
                    }
                }
            }
        }
    }
}
