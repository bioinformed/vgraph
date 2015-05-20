# Simple iterative graph traversal algorithms

from collections import deque


def dfs(graph, start):
    stack = [iter([start])]
    seen = set()
    while stack:
        node = next(stack[-1], None)
        if node is None:
            stack.pop()
        elif node not in seen:
            seen.add(node)
            stack.append(iter(graph[node]))
            yield node


def dfs_paths(graph, start):
    stack = [(iter([start]), [])]
    while stack:
        it, path = stack[-1]
        node = next(it, None)
        if node is None:
            stack.pop()
        elif node not in path:
            adj = graph[node]
            if adj:
                stack.append((iter(adj), path + [node]))
            else:
                yield path + [node]


def bfs(graph, start):
    queue = deque([start])
    seen = {start}
    while queue:
        node = queue.popleft()
        yield node
        queue += (n for n in graph[node] if n not in seen)
        seen.update(graph[node])


def bfs_paths(graph, start):
    queue = deque([[start]])
    while queue:
        path = queue.popleft()
        adj = [path + [node] for node in graph[path[-1]] if node not in path]
        if adj:
            queue += adj
        else:
            yield path
