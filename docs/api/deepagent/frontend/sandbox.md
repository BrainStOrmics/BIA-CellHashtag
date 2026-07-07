# Sandbox

[Frontend](/oss/python/deepagents/frontend/overview)[Patterns](/oss/python/deepagents/frontend/subagent-streaming)

# Sandbox
Copy page

Build an IDE-like UI for a coding agent backed by a sandbox environmentCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Coding agents need more than a chat window. They need a file browser, a code
viewer, and a diff panel, an IDE experience. This pattern connects a deep
agent to a [sandbox](/oss/python/deepagents/sandboxes) so it can read,
write, and execute code in an isolated environment, then exposes the sandbox
filesystem through a custom API server so the frontend can display files in
real time as the agent works.

## [​](#architecture)Architecture

The sandbox pattern has three layers:

{' ' * (self.list_depth - 1)}- 
**A deep agent with a sandbox backend:** The agent gets filesystem tools
(`read_file`, `write_file`, `edit_file`, `execute`) automatically from the
sandbox

{' ' * (self.list_depth - 1)}- 
**Custom API server** — A FastAPI app exposed via `langgraph.json`’s `http.app`
field, providing file browsing endpoints the frontend can call

{' ' * (self.list_depth - 1)}- 
**IDE frontend:** A three-panel layout (file tree, code/diff viewer, chat)
that syncs files in real time as the agent makes changes

## [​](#sandbox-lifecycle)Sandbox lifecycle

Before diving into the code, it’s important to understand how sandboxes are
scoped. The scoping strategy determines who shares a sandbox, how long it
lives, and how it’s resolved at runtime.

### [​](#thread-scoped-sandbox-recommended)Thread-scoped sandbox (recommended)

Each LangGraph thread gets its own sandbox. The sandbox ID is stored in the
thread’s metadata and resolved at runtime via `getConfig()`.
This is the recommended approach for most applications:

{' ' * (self.list_depth - 1)}- Conversations are isolated — file changes in one thread don’t affect another

{' ' * (self.list_depth - 1)}- Sandbox state persists across page reloads (same thread = same sandbox)

{' ' * (self.list_depth - 1)}- Cleanup is straightforward: when a thread is deleted, its sandbox can be too

### [​](#agent-scoped-sandbox)Agent-scoped sandbox

All threads under the same assistant share a single sandbox. Useful for
persistent project environments where you want changes to carry across
conversations:

```
from langgraph.config import get_config

def get_sandbox_backend_for_assistant():
 config = get_config()
 assistant_id = config.get("metadata", {}).get("assistant_id")
 return get_or_create_sandbox_for_assistant(assistant_id)

```

### [​](#user-scoped-sandbox)User-scoped sandbox

Each user gets their own sandbox across all threads. Requires custom
authentication and user identification:

```
from langgraph.config import get_config

def get_sandbox_backend_for_user():
 config = get_config()
 user_id = config.get("configurable", {}).get("user_id")
 return get_or_create_sandbox_for_user(user_id)

```

### [​](#session-scoped-sandbox-client-side)Session-scoped sandbox (client-side)

For simpler apps without LangGraph threads, the frontend can generate a
session ID and pass it directly. This approach doesn’t persist across
browser sessions and is best for demos or prototyping:

```
import uuid
import urllib.parse
import urllib.request

session_id = str(uuid.uuid4())
query = urllib.parse.urlencode({"sessionId": session_id})
urllib.request.urlopen(f"http://localhost:2024/api/sandbox/tree?{query}")

```

The rest of this guide uses **thread-scoped sandboxes** as the primary example.

## [​](#setting-up-the-agent)Setting up the agent

### [​](#choose-a-sandbox-provider)Choose a sandbox provider

Deep Agents supports multiple [sandbox providers](/oss/python/integrations/sandboxes). Any provider that implements the `SandboxBackendProtocol` works:

```
from deepagents import create_deep_agent
from deepagents.sandbox import LangSmithSandbox # or DaytonaSandbox, etc.

sandbox = LangSmithSandbox.create()
agent = create_deep_agent(model="google_genai:gemini-3.5-flash", backend=sandbox)

```

The agent automatically gets filesystem tools (`read_file`, `write_file`,
`edit_file`, `ls`, `glob`, `grep`) and an `execute` tool for running shell
commands. No tool configuration needed.

### [​](#resolve-a-sandbox-per-thread)Resolve a sandbox per thread

Instead of creating a sandbox at module level (which would be shared across
all threads and may expire), resolve the sandbox per-thread at runtime. The sandbox reads `thread_id` from the LangGraph config via `getConfig()`:

```
from deepagents import create_deep_agent
from deepagents.sandbox import LangSmithSandbox
from langgraph.config import get_config

def get_or_create_sandbox_for_thread(thread_id: str) -> LangSmithSandbox:
 # Look up or create sandbox based on thread_id
 ...

sandbox = LangSmithSandbox(
 resolve=lambda: get_or_create_sandbox_for_thread(
 get_config()["configurable"]["thread_id"]
 ),
)

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=sandbox,
)

```

### [​](#seed-the-sandbox)Seed the sandbox

Before the agent runs, populate the sandbox with your project files using
`uploadFiles`:
For **LangSmith** sandboxes, the container image and resource limits come from a
[sandbox snapshot](/langsmith/sandbox-snapshots). Pass `templateName` when creating
the sandbox (see `get_or_create_sandbox_for_thread` above). `upload_files` seeds or updatesproject files at runtime on top of that image.

```
const SEED_FILES: Record<string, string> = {
 "package.json": JSON.stringify({ name: "my-app", version: "1.0.0" }, null, 2),
 "src/index.js": 'console.log("Hello");',
};

const encoder = new TextEncoder();
await sandbox.uploadFiles(
 Object.entries(SEED_FILES).map(([path, content]) => [`/app/${path}`, encoder.encode(content)]),
);

```

Run `sandbox.execute("cd /app && npm install")` after uploading `package.json` to install
dependencies before the agent starts.

## [​](#adding-the-file-browsing-api)Adding the file browsing API

The agent can read and write files, but the frontend also needs direct access to
browse the sandbox filesystem. Add a custom [FastAPI](https://fastapi.tiangolo.com) API server
and expose it through the `http.app` field in `langgraph.json`.

### [​](#create-the-api-server)Create the API server

The sandbox API endpoints use the thread ID as a URL path parameter. This
ensures the frontend always accesses the correct sandbox for the current
conversation, using the same `get_or_create_sandbox_for_thread` function as the
agent’s backend:

```
# src/api/server.py
from fastapi import FastAPI, Query, Path
from utils import get_or_create_sandbox_for_thread

app = FastAPI()

@app.get("/api/sandbox/{thread_id}/tree")
async def list_tree(
 thread_id: str = Path(...),
 path: str = Query("/app"),
):
 sandbox = await get_or_create_sandbox_for_thread(thread_id)
 result = await sandbox.aexecute(
 f"find {path} -printf '%y\\t%s\\t%p\\n' 2>/dev/null | sort"
 )
 entries = []
 for line in result.output.strip().split("\n"):
 if not line:
 continue
 type_char, size_str, full_path = line.split("\t")
 entries.append({
 "name": full_path.split("/")[-1],
 "type": "directory" if type_char == "d" else "file",
 "path": full_path,
 "size": int(size_str),
 })
 return {"path": path, "entries": entries, "sandbox_id": sandbox.id}

@app.get("/api/sandbox/{thread_id}/file")
async def read_file(
 thread_id: str = Path(...),
 path: str = Query(...),
):
 sandbox = await get_or_create_sandbox_for_thread(thread_id)
 results = await sandbox.adownload_files([path])
 return {"path": path, "content": results[0].content.decode()}

```

Both the agent’s backend and the API server call the same
`get_or_create_sandbox_for_thread` function. This ensures they always resolveto the same sandbox for a given thread. The sandbox ID in thread metadata
is the single source of truth — no in-memory caches needed.

### [​](#configure-langgraph-json)Configure `langgraph.json`

Register both the agent graph and the API server. The `http.app` field tells
the LangGraph platform to serve your custom routes alongside the default ones:

```
{
 "graphs": {
 "coding_agent": "./src/agents/my_agent.py:agent"
 },
 "env": ".env",
 "http": {
 "app": "./src/api/server.py:app"
 }
}

```

Your custom routes are available at the same host as the LangGraph API. For
local development with `langgraph dev`, that’s `http://localhost:2024`.
Custom routes defined in `http.app` take priority over default LangGraph routes. This means you
can shadow built-in endpoints if needed, but be careful not to accidentally override routes like
`/threads` or `/runs`.

## [​](#building-the-frontend)Building the frontend

The frontend has three panels: a file tree sidebar, a code/diff viewer, and a
chat panel. It uses `useStream` for the agent conversation and the custom API
endpoints for file browsing.

### [​](#thread-creation)Thread creation

Create a LangGraph thread when the page loads and persist its ID in
`sessionStorage` so page reloads reconnect to the same sandbox:

```
const THREAD_KEY = "sandbox-thread-id";

function IDEPreview() {
 const [threadId, setThreadId] = useState<string | null>(
 () => sessionStorage.getItem(THREAD_KEY),
 );

 const updateThreadId = useCallback((id: string | null) => {
 setThreadId(id);
 if (id) sessionStorage.setItem(THREAD_KEY, id);
 else sessionStorage.removeItem(THREAD_KEY);
 }, []);

 const stream = useStream<typeof myAgent>({
 apiUrl: AGENT_URL,
 assistantId: "coding_agent",
 threadId,
 onThreadId: updateThreadId,
 });

 // Create thread on first mount
 useEffect(() => {
 if (threadId) return;
 stream.client.threads.create().then((t) => updateThreadId(t.thread_id));
 }, [stream.client, threadId, updateThreadId]);

 // Pass threadId to sandbox file hooks
 const { tree, files } = useSandboxFiles(threadId);
 // ...
}

```

The “new thread” button clears the stored ID so the next mount creates a
fresh thread (and sandbox):

```
function handleNewThread() {
 stream.switchThread(null);
 updateThreadId(null);
}

```

### [​](#file-state-management)File state management

Track two snapshots of the sandbox filesystem: the original state (before the
agent runs) and the current state (updated in real time). The thread ID is
included in the API URL so requests always hit the correct sandbox:

```
const AGENT_URL = "http://localhost:2024";

async function fetchTree(threadId: string): Promise<FileEntry[]> {
 const res = await fetch(
 `${AGENT_URL}/api/sandbox/${encodeURIComponent(threadId)}/tree?filePath=/app`,
 );
 const data = await res.json();
 return data.entries.filter((e: FileEntry) => !e.path.includes("node_modules"));
}

async function fetchFile(threadId: string, path: string): Promise<string | null> {
 const res = await fetch(
 `${AGENT_URL}/api/sandbox/${encodeURIComponent(threadId)}/file?filePath=${encodeURIComponent(path)}`,
 );
 const data = await res.json();
 return data.content ?? null;
}

```

### [​](#real-time-file-sync)Real-time file sync

The key to the IDE experience is updating files **as the agent works**, not
after it finishes. Watch the stream’s messages for `ToolMessage` instances
from file-mutating tools. When a `write_file` or `edit_file` tool call
completes, refresh that specific file. When `execute` completes, refresh
everything (since a shell command could modify any file):
ReactVueSvelteAngular
```
import { useStream } from "@langchain/react";
import { ToolMessage, AIMessage } from "langchain";

const FILE_MUTATING_TOOLS = new Set(["write_file", "edit_file", "execute"]);

export function IDEPreview() {
 const stream = useStream<typeof myAgent>({
 apiUrl: AGENT_URL,
 assistantId: "coding_agent",
 });

 const processedIds = useRef(new Set<string>());

 useEffect(() => {
 // Build a map of file-mutating tool calls from AI messages
 const toolCallMap = new Map();
 for (const msg of stream.messages) {
 if (!AIMessage.isInstance(msg)) continue;
 for (const tc of msg.tool_calls ?? []) {
 if (tc.id && FILE_MUTATING_TOOLS.has(tc.name)) {
 toolCallMap.set(tc.id, { name: tc.name, args: tc.args });
 }
 }
 }

 // When a ToolMessage appears for a file-mutating tool, refresh
 for (const msg of stream.messages) {
 if (!ToolMessage.isInstance(msg)) continue;
 const id = msg.id ?? msg.tool_call_id;
 if (!id || processedIds.current.has(id)) continue;

 const call = toolCallMap.get(msg.tool_call_id);
 if (!call) continue;
 processedIds.current.add(id);

 if (call.name === "write_file" || call.name === "edit_file") {
 refreshSingleFile(call.args.path);
 } else if (call.name === "execute") {
 refreshAllFiles();
 }
 }
 }, [stream.messages]);
}

```

### [​](#detecting-changed-files)Detecting changed files

Before each agent run, snapshot the current file contents. After files refresh,
compare against the snapshot to identify which files changed:

```
function detectChanges(current: FileSnapshot, original: FileSnapshot): Set<string> {
 const changed = new Set<string>();
 for (const [path, content] of Object.entries(current)) {
 if (original[path] !== content) changed.add(path);
 }
 for (const path of Object.keys(original)) {
 if (!(path in current)) changed.add(path);
 }
 return changed;
}

```

When a user selects a changed file, default to the diff view so they
immediately see what the agent modified.

### [​](#displaying-diffs)Displaying diffs

Use a framework-appropriate diff library to render unified diffs:

| 
 | Framework | Library | Component
 | React | [`@pierre/diffs`](https://diffs.com) | `<FileDiff>` with `parseDiffFromFile`
 | Vue | [`@git-diff-view/vue`](https://github.com/MrWangJustToDo/git-diff-view) | `<DiffView>` with `generateDiffFile` from `@git-diff-view/file`
 | Svelte | [`@git-diff-view/svelte`](https://github.com/MrWangJustToDo/git-diff-view) | `<DiffView>` with `generateDiffFile` from `@git-diff-view/file`
 | Angular | [`ngx-diff`](https://github.com/rars/ngx-diff) | `<ngx-unified-diff>` with `[before]` and `[after]`
Example with `@pierre/diffs` (React):

```
import { FileDiff } from "@pierre/diffs/react";
import { parseDiffFromFile } from "@pierre/diffs";

function DiffPanel({ original, current, fileName }) {
 const diff = parseDiffFromFile(
 { name: fileName, contents: original },
 { name: fileName, contents: current },
 );

 return (
 <FileDiff
 fileDiff={diff}
 options={{ theme: "github-dark", diffStyle: "unified", diffIndicators: "bars" }}
 />
 );
}

```

### [​](#changed-files-summary)Changed files summary

Show a summary of all modified files with line-level addition/deletion counts.
This gives users a quick overview of the agent’s impact — similar to a `git status`:

```
function ChangedFilesSummary({ changedFiles, files, originalFiles, onSelect }) {
 const stats = [...changedFiles].map((path) => {
 const oldLines = (originalFiles[path] ?? "").split("\n");
 const newLines = (files[path] ?? "").split("\n");
 // Compute additions/deletions by comparing lines
 return { path, additions, deletions };
 });

 return (
 <div>
 <h3>{stats.length} Files Changed</h3>
 {stats.map((file) => (
 <button key={file.path} onClick={() => onSelect(file.path)}>
 {file.path}
 <span className="text-green-400">+{file.additions}</span>
 <span className="text-red-400">-{file.deletions}</span>
 </button>
 ))}
 </div>
 );
}

```

## [​](#the-three-panel-layout)The three-panel layout

The IDE layout arranges three panels side by side:

| 
 | Panel | Width | Purpose
 | File tree | Fixed (208px) | Browse sandbox files, see change indicators
 | Code / Diff | Flexible | View file content or unified diff
 | Chat | Fixed (320px) | Interact with the agent

```
<div className="flex h-screen">
 <div className="w-52 shrink-0">
 <FileTree />
 <ChangedFilesSummary />
 </div>

 <CodePanel /* flex-1 */ />

 <div className="w-80 shrink-0">
 <ChatPanel />
 </div>
</div>

```

The file tree shows VS Code-style icons (using
[`@iconify-json/vscode-icons`](https://www.npmjs.com/package/@iconify-json/vscode-icons))
and amber dots on modified files. Selecting a modified file automatically
switches to the diff tab.

## [​](#use-cases)Use cases

A sandbox is the right choice when:

{' ' * (self.list_depth - 1)}- **Coding agents** that create, modify, and run code need a visual interface
beyond chat

{' ' * (self.list_depth - 1)}- **Code review workflows** where the agent suggests changes and the user
reviews diffs before accepting

{' ' * (self.list_depth - 1)}- **Tutorial or learning apps** where an AI assistant helps users build a
project step by step, showing changes in context

{' ' * (self.list_depth - 1)}- **Prototyping tools** where users describe features in natural language and
watch the agent implement them in real time

## [​](#best-practices)Best practices

{' ' * (self.list_depth - 1)}- **Use thread-scoped sandboxes** for production apps. Store the sandbox ID in
thread metadata and resolve it via `getConfig()` at runtime. This avoids
module-level state and keeps sandboxes isolated per conversation.

{' ' * (self.list_depth - 1)}- **Share `getOrCreateSandboxForThread`** between the agent backend and the
API server. Both should resolve the sandbox the same way — via thread
metadata — so there’s a single source of truth with no in-memory caches.

{' ' * (self.list_depth - 1)}- **Persist `threadId` in `sessionStorage`** so page reloads reconnect to the
same thread and sandbox instead of creating new ones.

{' ' * (self.list_depth - 1)}- **Sync files on every relevant tool call**, not just when the run finishes. This
makes the IDE feel live. Watch for `write_file`, `edit_file`, and `execute`
tool messages and refresh immediately.

{' ' * (self.list_depth - 1)}- **Default to diff view for changed files**. When a user clicks a file that
was modified by the agent, show the diff first — that’s what they care about.

{' ' * (self.list_depth - 1)}- **Show compact tool results for read-only operations**. Instead of dumping
the full output of `read_file` in the chat, show a one-liner like
`Read router.js L1-42`. Reserve the full output display for mutating tools.

{' ' * (self.list_depth - 1)}- **Seed the sandbox with a real project**. Starting from an empty sandbox is
disorienting. Upload a working starter project so users (and the agent) have
context immediately.

{' ' * (self.list_depth - 1)}- **Filter `node_modules` from the file tree**. Nobody wants to browse
thousands of dependency files. Filter them out when fetching the tree.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/deepagents/frontend/sandbox.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Todo listPrevious](/oss/python/deepagents/frontend/todo-list)[Agent Client Protocol (ACP)Next](/oss/python/deepagents/acp)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
