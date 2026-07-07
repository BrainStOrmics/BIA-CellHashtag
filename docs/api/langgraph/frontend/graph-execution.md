# Graph execution

[Frontend](/oss/python/langgraph/frontend/overview)

# Graph execution
Copy page

Visualize multi-step graph pipelines with per-node status and streaming contentCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.LangGraph agents aren’t black boxes. Every graph is composed of **named nodes**
that execute in sequence or in parallel: classify, research, analyze,
synthesize. Graph execution cards make this pipeline visible by rendering a card
for each node, showing its status, streaming its content in real time, and
tracking completion across the entire workflow. Users see exactly what the agent
is doing, which step it’s on, and what each step produced.

## [​](#how-graph-nodes-map-to-ui-cards)How graph nodes map to UI cards

A LangGraph graph defines a series of nodes, each responsible for a specific
task. For example, a research pipeline might have:

{' ' * (self.list_depth - 1)}- **Classify**: categorize the user’s query

{' ' * (self.list_depth - 1)}- **Research**: gather relevant information

{' ' * (self.list_depth - 1)}- **Analyze**: draw conclusions from the research

{' ' * (self.list_depth - 1)}- **Synthesize**: produce a final, polished response

Each node writes its output to a specific key in the graph’s state. By mapping
these node names and state keys to UI components, you can create a visual
representation of the entire pipeline.

```
const PIPELINE_NODES = [
 { name: "classify", stateKey: "classification", label: "Classify" },
 { name: "do_research", stateKey: "research", label: "Research" },
 { name: "analyze", stateKey: "analysis", label: "Analyze" },
 { name: "synthesize", stateKey: "synthesis", label: "Synthesize" },
];

const PIPELINE_NODE_NAMES = new Set(PIPELINE_NODES.map((n) => n.name));

```

## [​](#setting-up-usestream)Setting up useStream

Wire up `useStream` as usual. The key properties you’ll use are `messages`
(for streaming content routing), `values` (for completed node outputs), and
`getMessagesMetadata` (for identifying which node produced each token).
Define a TypeScript interface matching your agent’s state schema and pass it as a type parameter to `useStream` for type-safe access to state values, including custom state keys for each pipeline node. In the examples below, replace `typeof myAgent` with your interface name:

```
import type { BaseMessage } from "@langchain/core/messages";

interface AgentState {
 messages: BaseMessage[];
 classification: string;
 research: string;
 analysis: string;
 synthesis: string;
}

```

ReactVueSvelteAngular
```
import { useStream } from "@langchain/react";

const AGENT_URL = "http://localhost:2024";

export function PipelineChat() {
 const stream = useStream<typeof myAgent>({
 apiUrl: AGENT_URL,
 assistantId: "graph_execution_cards",
 });

 return (
 <div>
 <PipelineProgress nodes={PIPELINE_NODES} values={stream.values} />
 <NodeCardList
 nodes={PIPELINE_NODES}
 messages={stream.messages}
 values={stream.values}
 getMetadata={stream.getMessagesMetadata}
 />
 </div>
 );
}

```

## [​](#routing-streaming-tokens-to-nodes)Routing streaming tokens to nodes

As the agent streams, each message is annotated with metadata that identifies
which graph node produced it. Use `getMessagesMetadata` to extract the
`langgraph_node` value and route tokens to the correct card:

```
function getStreamingContent(
 messages: BaseMessage[],
 getMetadata: (msg: BaseMessage) => MessageMetadata | undefined
): Record<string, string> {
 const content: Record<string, string> = {};

 for (const message of messages) {
 if (message.type !== "ai") continue;

 const metadata = getMetadata(message);
 const node = metadata?.streamMetadata?.langgraph_node;

 if (node && PIPELINE_NODE_NAMES.has(node)) {
 content[node] = typeof message.content === "string"
 ? message.content
 : "";
 }
 }

 return content;
}

```

This gives you a map from node name to its current streaming content. As tokens
arrive, the corresponding card updates in real time.
The `streamMetadata.langgraph_node` field is set automatically by LangGraph.
You don’t need any special configuration on the backend. Just stream messages
as usual, and the metadata is included.

## [​](#determining-node-status)Determining node status

Each node can be in one of four states: not started, streaming, complete, or
idle. You derive the status from two sources: the streaming content map (for
active nodes) and `stream.values` (for completed nodes):

```
type NodeStatus = "idle" | "streaming" | "complete";

function getNodeStatus(
 node: { name: string; stateKey: string },
 streamingContent: Record<string, string>,
 values: Record<string, unknown>
): NodeStatus {
 if (values?.[node.stateKey]) return "complete";
 if (streamingContent[node.name]) return "streaming";
 return "idle";
}

```

## [​](#building-the-pipeline-progress-bar)Building the pipeline progress bar

A horizontal progress bar at the top gives users a bird’s-eye view of the
entire pipeline. Each step is a labeled segment that fills in as nodes complete:

```
function PipelineProgress({
 nodes,
 values,
 streamingContent,
}: {
 nodes: typeof PIPELINE_NODES;
 values: Record<string, unknown>;
 streamingContent: Record<string, string>;
}) {
 return (
 <div className="flex items-center gap-1">
 {nodes.map((node, i) => {
 const status = getNodeStatus(node, streamingContent, values);
 const colors = {
 idle: "bg-gray-200 text-gray-500",
 streaming: "bg-blue-400 text-white animate-pulse",
 complete: "bg-green-500 text-white",
 };

 return (
 <div key={node.name} className="flex items-center">
 <div
 className={`rounded-full px-3 py-1 text-xs font-medium ${colors[status]}`}
 >
 {node.label}
 </div>
 {i < nodes.length - 1 && (
 <div
 className={`mx-1 h-0.5 w-6 ${
 status === "complete" ? "bg-green-500" : "bg-gray-200"
 }`}
 />
 )}
 </div>
 );
 })}
 </div>
 );
}

```

## [​](#building-collapsible-nodecard-components)Building collapsible NodeCard components

Each node gets its own card that shows the status badge, content (streaming or
final), and a collapsible body for long outputs:

```
function NodeCard({
 node,
 status,
 streamingContent,
 completedContent,
}: {
 node: { name: string; stateKey: string; label: string };
 status: NodeStatus;
 streamingContent: string | undefined;
 completedContent: unknown;
}) {
 const [collapsed, setCollapsed] = useState(false);

 const displayContent =
 status === "complete"
 ? formatContent(completedContent)
 : streamingContent ?? "";

 const statusBadge = {
 idle: { text: "Waiting", className: "bg-gray-100 text-gray-600" },
 streaming: {
 text: "Running",
 className: "bg-blue-100 text-blue-700 animate-pulse",
 },
 complete: { text: "Done", className: "bg-green-100 text-green-700" },
 };

 const badge = statusBadge[status];

 return (
 <div className="rounded-lg border bg-white shadow-sm">
 <button
 onClick={() => setCollapsed(!collapsed)}
 className="flex w-full items-center justify-between p-4"
 >
 <div className="flex items-center gap-3">
 <h3 className="font-semibold">{node.label}</h3>
 <span
 className={`rounded-full px-2 py-0.5 text-xs font-medium ${badge.className}`}
 >
 {badge.text}
 </span>
 </div>
 <ChevronIcon direction={collapsed ? "down" : "up"} />
 </button>

 {!collapsed && displayContent && (
 <div className="border-t px-4 py-3">
 <div className="prose prose-sm max-w-none">
 {displayContent}
 {status === "streaming" && (
 <span className="inline-block h-4 w-1 animate-pulse bg-blue-500" />
 )}
 </div>
 </div>
 )}
 </div>
 );
}

function formatContent(value: unknown): string {
 if (typeof value === "string") return value;
 if (value == null) return "";
 return JSON.stringify(value, null, 2);
}

```

## [​](#streaming-vs-completed-content)Streaming vs. completed content

There are two sources of content for each node, and picking the right one
matters for a smooth UX:

| 
 | Source | When to use
 | `streamingContent[node.name]` | While the node is actively streaming, this contains tokens as they arrive
 | `stream.values[node.stateKey]` | After the node completes, this contains the final, committed output
The pattern is: show streaming content for live updates, fall back to the
committed state value once the node is done.

```
for (const node of PIPELINE_NODES) {
 const status = getNodeStatus(node, streamingContent, stream.values);

 const content =
 status === "streaming"
 ? streamingContent[node.name]
 : stream.values?.[node.stateKey];
}

```

Streaming content may include partial tokens or markdown that hasn’t been
fully formed yet. If you render markdown, make sure your renderer handles
incomplete syntax gracefully (e.g., an unclosed bold marker `**`).

## [​](#putting-it-all-together)Putting it all together

Here’s the full card list that combines routing, status detection, and card
rendering:

```
function NodeCardList({
 nodes,
 messages,
 values,
 getMetadata,
}: {
 nodes: typeof PIPELINE_NODES;
 messages: BaseMessage[];
 values: Record<string, unknown>;
 getMetadata: (msg: BaseMessage) => MessageMetadata | undefined;
}) {
 const streamingContent = getStreamingContent(messages, getMetadata);

 return (
 <div className="space-y-3">
 {nodes.map((node) => {
 const status = getNodeStatus(node, streamingContent, values);
 return (
 <NodeCard
 key={node.name}
 node={node}
 status={status}
 streamingContent={streamingContent[node.name]}
 completedContent={values?.[node.stateKey]}
 />
 );
 })}
 </div>
 );
}

```

## [​](#use-cases)Use cases

Graph execution cards work well for any multi-step pipeline where visibility
matters:

{' ' * (self.list_depth - 1)}- **Research pipelines**: classify → gather sources → analyze → synthesize a
report

{' ' * (self.list_depth - 1)}- **Content generation**: outline → draft → fact-check → edit → publish

{' ' * (self.list_depth - 1)}- **Data processing**: ingest → validate → transform → aggregate → export

{' ' * (self.list_depth - 1)}- **Code generation**: understand requirements → plan architecture → write
code → review → test

{' ' * (self.list_depth - 1)}- **Decision workflows**: gather context → evaluate options → score
alternatives → recommend

## [​](#handling-dynamic-pipelines)Handling dynamic pipelines

Not all graphs have a fixed set of nodes. Some pipelines add or skip nodes
based on the input. Handle this by checking which state keys actually have
values:

```
const activeNodes = PIPELINE_NODES.filter(
 (node) =>
 streamingContent[node.name] ||
 values?.[node.stateKey] ||
 node.name === currentNode
);

```

This ensures your UI only shows cards for nodes that are relevant to the
current execution, avoiding empty placeholder cards.
If your graph has conditional branching (e.g., skip “Research” for simple
factual queries), the skipped nodes will never appear in the streaming content
or state values. Your pipeline progress bar should reflect this by dimming or
hiding skipped steps.

## [​](#best-practices)Best practices

{' ' * (self.list_depth - 1)}- **Define nodes declaratively**. Keep your `PIPELINE_NODES` array as a single
source of truth that maps node names, state keys, and display labels.

{' ' * (self.list_depth - 1)}- **Prefer streaming content for active nodes**. It gives users immediate
feedback. Only fall back to committed state values after the node completes.

{' ' * (self.list_depth - 1)}- **Auto-collapse completed nodes**. In long pipelines, auto-collapse finished
cards so users can focus on the currently active step.

{' ' * (self.list_depth - 1)}- **Show estimated timing**. If you have historical data on how long each node
takes, display a time estimate to set user expectations.

{' ' * (self.list_depth - 1)}- **Add a global progress indicator**. Complement per-node cards with an
overall progress bar (e.g., “Step 2 of 4”) at the top of the pipeline view.

{' ' * (self.list_depth - 1)}- **Handle errors per node**. If a node fails, show the error in its card
without collapsing the entire pipeline. Other nodes may still complete
successfully.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langgraph/frontend/graph-execution.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[OverviewPrevious](/oss/python/langgraph/frontend/overview)[Choosing between Graph and Functional APIsNext](/oss/python/langgraph/choosing-apis)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
