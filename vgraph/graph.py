# Copyright 2015 Kevin B Jacobs
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may
# not use this file except in compliance with the License.  You may obtain
# a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  See the
# License for the specific language governing permissions and limitations
# under the License.

"""Simple iterative graph traversal algorithms."""


from collections import deque


def dfs(graph, start):
    """Depth first search."""
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
    """Depth first search yielding the path taken."""
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
    """Breadth first search."""
    queue = deque([start])
    seen = {start}
    while queue:
        node = queue.popleft()
        yield node
        queue += (n for n in graph[node] if n not in seen)
        seen.update(graph[node])


def bfs_paths(graph, start):
    """Breadth first search yielding the path taken."""
    queue = deque([[start]])
    while queue:
        path = queue.popleft()
        adj = [path + [node] for node in graph[path[-1]] if node not in path]
        if adj:
            queue += adj
        else:
            yield path
