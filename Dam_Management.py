import heapq
import tkinter as tk
import copy
import matplotlib.pyplot as plt
import networkx as nx
import time


class GraphVisualization:
    def __init__(self):
        self.visual = []
        self.edge_labels = {}
        
    def addEdge(self, a, b, i):
        temp = [a, b]
        self.visual.append(temp)
        self.edge_labels[(a, b)] = i
        
    def visualize(self):
        G = nx.DiGraph()  # Use nx.DiGraph() to create a directed graph
        G.add_edges_from(self.visual)
        
        pos = nx.spring_layout(G)
        nx.draw_networkx(G, pos, arrows=True)  # Set arrows=True to show directed edges
        nx.draw_networkx_edge_labels(G, pos, edge_labels=self.edge_labels, font_color='red')
        
        plt.axis('off')
        plt.show()


G = GraphVisualization()
def generate_graph(adjacency_list):
    for idx, adj in enumerate(adjacency_list):
        for idx2, i in enumerate(adj):
            if i != 0:
                if adjacency_list[idx2][idx]!=0:
                    G.addEdge(idx,idx2,[adjacency_list[idx2][idx],i])
                else:
                    G.addEdge(idx, idx2, i)
    G.visualize()

paths_list=[]
def highlight(lis):
    global paths_list
    paths_list.append(lis)
    print(paths_list)

class WaterResourceManagementSystem:
    def __init__(self, graph, total_volume):
        self.graph = graph
        self.num_d = len(graph)
        self.total_volume = total_volume

    def dijkstra(self, source, sink):
        dist = [float("inf")] * self.num_d
        dist[source] = 0
        pq = [(0, source)]
        while pq:
            (d, u) = heapq.heappop(pq)
            if d > dist[u]:
                continue
            for v, w in enumerate(self.graph[u]):
                if w > 0 and dist[u] + w < dist[v]:
                    dist[v] = dist[u] + w
                    heapq.heappush(pq, (dist[v], v))
        return dist[sink]

    def bfs(self, source, sink, parent):
        visited = [False] * self.num_d
        queue = []
        queue.append(source)
        visited[source] = True
        while queue:
            u = queue.pop(0)
            for ind, val in enumerate(self.graph[u]):
                if visited[ind] == False and val > 0:
                    queue.append(ind)
                    visited[ind] = True
                    parent[ind] = u
        return True if visited[sink] else False

    def ford_fulkerson(self, source, sink):
        parent = [-1] * self.num_d
        max_flow = 0
        current_path = []  # Initialize the current_path variable
        while self.bfs(source, sink, parent):
            sss=[]
            path_flow = float("Inf")
            s = sink
            while s != source:
                path_flow = min(path_flow, self.graph[parent[s]][s])
                s = parent[s]
                sss.append(s)
            max_flow += path_flow
            v = sink
            while v != source:
                u = parent[v]
                self.graph[u][v] -= path_flow
                self.graph[v][u] += path_flow
                v = parent[v]
            rev=sss[::-1]
            rev.append(len(self.graph)-1)
            highlight(rev)
        return max_flow
    
    def bfs_dinic(self, source, sink, level):
            queue = []
            queue.append(source)
            level[source] = 1
            while queue:
                u = queue.pop(0)
                for v, w in enumerate(self.graph[u]):
                    if level[v] == 0 and w > 0:
                        queue.append(v)
                        level[v] = level[u] + 1
            return True if level[sink] > 0 else False

    def send_flow_dinic(self, u, flow, sink, level, blocking_flows):
        if u == sink:
            return flow
        while blocking_flows[u] < len(self.graph[u]):
            v = blocking_flows[u]
            if self.graph[u][v] > 0 and level[v] == level[u] + 1:
                min_flow = min(flow, self.graph[u][v])
                bottleneck_flow = self.send_flow_dinic(v, min_flow, sink, level, blocking_flows)
                if bottleneck_flow > 0:
                    self.graph[u][v] -= bottleneck_flow
                    self.graph[v][u] += bottleneck_flow
                    return bottleneck_flow
            blocking_flows[u] += 1
        return 0

    def dinic(self, source, sink):
        max_flow = 0
        level = [0] * self.num_d
        while self.bfs_dinic(source, sink, level):
            blocking_flows = [0] * self.num_d
            while True:
                flow = self.send_flow_dinic(source, float("inf"), sink, level, blocking_flows)
                if flow == 0:
                    break
                max_flow += flow
        return max_flow
    
    def findMaxVol(self,source,sink):
        s=0
        graph_copy = copy.deepcopy(self.graph)
        q=[]
        q.append(graph_copy[source])
        ret=[]
        visited=[]
        visited.append(source)
        while q:
            node=q.pop(0)
            for i in range(len(node)):
                if node[i]!=0:
                    s+=node[i]
                    if i not in visited:
                        visited.append(i)
                        q.append(graph_copy[i])
                        ret.append(i)
        return s,ret



    def find_overflow(self, source, sink, vol):
        max_vol,lis = self.findMaxVol(source,sink)
        if max_vol<vol:
            return (True,lis,max_vol)
        else:
            return (False,lis,max_vol)


class WaterResourceManagementSystemGUI:
    def __init__(self, master):
        self.master = master
        self.master.title("Water Resource Management System")
        self.master.geometry("500x700")

        # Set font family and size
        self.font_family = "Arial"
        self.font_size = 12

        # Set colors
        self.bg_color = "#f2cc8a"
        self.button_color = "#032d69"
        self.button_text_color = "white"

        # Set padding and spacing
        self.padding = 7
        self.spacing = 7

        # Set the background color of the window
        self.master.config(bg=self.bg_color)

        # Create input fields for the graph
        self.graph_label = tk.Label(self.master, text="Enter the graph as a list of lists:", font=(self.font_family, self.font_size), bg=self.bg_color)
        self.graph_label.pack(pady=self.padding, padx=self.padding)
        self.graph_entry = tk.Entry(self.master, width=50, font=(self.font_family, self.font_size))
        self.graph_entry.pack(pady=self.spacing, padx=self.padding)

        # Create input field for the total water volume
        self.volume_label = tk.Label(self.master, text="Enter the total water volume:", font=(self.font_family, self.font_size), bg=self.bg_color)
        self.volume_label.pack(pady=self.padding, padx=self.padding)
        self.volume_entry = tk.Entry(self.master, width=50, font=(self.font_family, self.font_size))
        self.volume_entry.pack(pady=self.spacing, padx=self.padding)

        # Create input fields for the source and sink
        self.source_label = tk.Label(self.master, text="Enter the source:", font=(self.font_family, self.font_size), bg=self.bg_color)
        self.source_label.pack(pady=self.padding, padx=self.padding)
        self.source_entry = tk.Entry(self.master, width=50, font=(self.font_family, self.font_size))
        self.source_entry.pack(pady=self.spacing, padx=self.padding)

        self.sink_label = tk.Label(self.master, text="Enter the sink:", font=(self.font_family, self.font_size), bg=self.bg_color)
        self.sink_label.pack(pady=self.padding, padx=self.padding)
        self.sink_entry = tk.Entry(self.master, width=50, font=(self.font_family, self.font_size))
        self.sink_entry.pack(pady=self.spacing, padx=self.padding)

        # Create a button to run the algorithm
        self.run_button = tk.Button(self.master, text="Run", command=self.run_algorithm, bg=self.button_color, fg=self.button_text_color, padx=10, pady=5, font=(self.font_family, self.font_size))
        self.run_button.pack(pady=self.padding, padx=self.padding)

        self.run_button = tk.Button(self.master, text="Show", command=self.show, bg=self.button_color, fg=self.button_text_color, padx=10, pady=5, font=(self.font_family, self.font_size))
        self.run_button.pack(pady=self.padding, padx=self.padding)

        self.run_button = tk.Button(self.master, text="Paths", command=self.path, bg=self.button_color, fg=self.button_text_color, padx=10, pady=5, font=(self.font_family, self.font_size))
        self.run_button.pack(pady=self.padding, padx=self.padding)

        # Create a label to display the output of the algorithm
        self.output_label = tk.Label(self.master, text="", font=(self.font_family, self.font_size), bg=self.bg_color)
        self.output_label.pack(pady=self.padding, padx=self.padding)

        # Create a button to find the shortest path
        self.shortest_path_button = tk.Button(self.master, text="Find Shortest Path", command=self.find_shortest_path, bg=self.button_color, fg=self.button_text_color, padx=10, pady=5, font=(self.font_family, self.font_size))
        self.shortest_path_button.pack(pady=self.padding, padx=self.padding)

        # Create a label to display the output of the shortest path algorithm
        self.shortest_path_label = tk.Label(self.master, text="", font=(self.font_family, self.font_size), bg=self.bg_color)
        self.shortest_path_label.pack(pady=self.padding, padx=self.padding)

        # Create a button to find overflow and max volume
        self.overflow_button = tk.Button(self.master, text="Find Overflow and Max Volume", command=self.find_overflow, bg=self.button_color, fg=self.button_text_color, padx=10, pady=5, font=(self.font_family, self.font_size))
        self.overflow_button.pack(pady=self.padding, padx=self.padding)

        # Create a label to display the output of the overflow algorithm
        self.overflow_label = tk.Label(self.master, text="", font=(self.font_family, self.font_size), bg=self.bg_color)
        self.overflow_label.pack(pady=self.padding, padx=self.padding)

    def path(self):
        global paths_list
        self.output_label.config(text=f"Paths: {paths_list}")

    def show(self):
        graph_str = self.graph_entry.get()
        graph = eval(graph_str)
        generate_graph(graph)

    def run_algorithm(self):
        # Parse the input fields
        graph_str = self.graph_entry.get()
        volume_str = self.volume_entry.get()
        source_str = self.source_entry.get()
        sink_str = self.sink_entry.get()

        # Convert the graph input to a list of lists
        graph = eval(graph_str)
        total_volume = int(volume_str)
        source = int(source_str)
        sink = int(sink_str)

        # Create an instance of the WaterResourceManagementSystem class
        wrms = WaterResourceManagementSystem(graph, total_volume)

        # Run the Ford-Fulkerson algorithm and measure the time taken
        start_time_ffa = time.perf_counter()
        max_flow_ffa = wrms.ford_fulkerson(source, sink)
        end_time_ffa = time.perf_counter()
        time_taken_ffa = (end_time_ffa - start_time_ffa) * 1000 

        # Run the Dinic's algorithm and measure the time taken
        start_time_dinic = time.perf_counter()
        max_flow_dinic = wrms.dinic(source, sink)
        end_time_dinic = time.perf_counter()
        time_taken_dinic = (end_time_dinic - start_time_dinic) * 1000

        # Display the output
        self.output_label.config(
            text=f"Max Flow (FFA): {max_flow_ffa}\nTime taken (FFA): {time_taken_ffa} seconds"
            f"\n\nMax Flow (Dinic's): {max_flow_ffa}\nTime taken (Dinic's): {time_taken_dinic} seconds"
        )


    def find_shortest_path(self):
        # Parse the input fields
        graph_str = self.graph_entry.get()
        source_str = self.source_entry.get()
        sink_str = self.sink_entry.get()

        # Convert the graph input to a list of lists
        graph = eval(graph_str)
        source = int(source_str)
        sink = int(sink_str)

        # Create an instance of the WaterResourceManagementSystem class
        wrms = WaterResourceManagementSystem(graph, 0)

        # Run Dijkstra's algorithm
        shortest_path = wrms.dijkstra(source, sink)

        # Display the output
        self.shortest_path_label.config(text=f"Shortest Path: {shortest_path}")

    def find_overflow(self):
    # Parse the input fields
        graph_str = self.graph_entry.get()
        volume_str = self.volume_entry.get()
        source_str = self.source_entry.get()
        sink_str = self.sink_entry.get()

        # Convert the graph input to a list of lists
        graph = eval(graph_str)
        total_volume = int(volume_str)
        source = int(source_str)
        sink = int(sink_str)

        # Create an instance of the WaterResourceManagementSystem class
        wrms = WaterResourceManagementSystem(graph, total_volume)

        # Run the overflow algorithm
        overflow, overflow_dams, max_volume = wrms.find_overflow(source, sink, total_volume)

        # Display the output
        if overflow:
            overflow_dams_str = ", ".join(str(dam) for dam in overflow_dams)
            self.overflow_label.config(text=f"Overflow \n Dams: {overflow_dams_str}\nMax Volume: {max_volume}")
        else:
            self.overflow_label.config(text=f"No overflow \nMax Volume: {max_volume}")


# Create the Tkinter application
root = tk.Tk()

# Create an instance of the WaterResourceManagementSystemGUI class
wrms_gui = WaterResourceManagementSystemGUI(root)

# Run the Tkinter event loop
root.mainloop()