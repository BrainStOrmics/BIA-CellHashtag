# Subagent streaming

[Frontend](/oss/python/deepagents/frontend/overview)[Patterns](/oss/python/deepagents/frontend/subagent-streaming)

# Subagent streaming
Copy page

Display specialist subagents with streaming content, progress tracking, and collapsible cardsCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.When a coordinator agent spawns specialist subagents (a researcher, an
analyst, a writer), you need to render the orchestrator’s messages separately
from each subagent’s streaming output. Set `filterSubagentMessages: true` in
`useStream` to cleanly split these two streams, then use
`getSubagentsByMessage` to attach each subagent’s progress card to the
coordinator message that triggered it.

## [​](#why-filter-subagent-messages)Why filter subagent messages

Without filtering, every token produced by every subagent appears interleaved
in the coordinator’s message stream, making it unreadable. With
`filterSubagentMessages: true`:

{' ' * (self.list_depth - 1)}- `stream.messages` contains only the coordinator’s messages

{' ' * (self.list_depth - 1)}- Each subagent’s content is accessible through `stream.subagents` and
`stream.getSubagentsByMessage`

{' ' * (self.list_depth - 1)}- The UI stays clean: the coordinator’s reasoning is separate from the
specialists’ work

This separation lets you render the orchestrator’s messages in one place and
attach each subagent’s progress card exactly where it belongs: beneath the
coordinator message that spawned it.

## [​](#setting-up-usestream)Setting up useStream

Always set `filterSubagentMessages: true`. This removes subagent tokens from
the main message stream so you can render the coordinator’s messages and
subagent output independently.
Define a TypeScript interface matching your agent’s state schema and pass it as a type parameter to `useStream` for type-safe access to state values. In the examples below, replace `typeof myAgent` with your interface name:

```
import type { BaseMessage } from "@langchain/core/messages";

interface AgentState {
 messages: BaseMessage[];
}

```

ReactVueSvelteAngular
```
import { useStream } from "@langchain/react";

const AGENT_URL = "http://localhost:2024";

export function DeepAgentChat() {
 const stream = useStream<typeof myAgent>({
 apiUrl: AGENT_URL,
 assistantId: "deep_agent_subagent_cards",
 filterSubagentMessages: true,
 });

 return (
 <div>
 {stream.messages.map((msg) => (
 <MessageWithSubagents
 key={msg.id}
 message={msg}
 subagents={stream.getSubagentsByMessage(msg.id)}
 />
 ))}
 </div>
 );
}

```

## [​](#submitting-with-subgraph-streaming)Submitting with subgraph streaming

When submitting a message, enable subgraph streaming and set an appropriate
recursion limit. Deep agent workflows often involve multiple layers of nested
subagraphs, so a higher recursion limit prevents premature termination:

```
stream.submit(
 { messages: [{ type: "human", content: text }] },
 { streamSubgraphs: true }
);

```

Deep Agents sets a default recursion limit of 10,000, which is sufficient for
most multi-expert setups. You can override this via `config.recursion_limit` if
needed.

## [​](#the-subagentstreaminterface)The SubagentStreamInterface

Each subagent exposes a `SubagentStreamInterface` with metadata about the
subagent’s task, status, and timing:

```
interface SubagentStreamInterface {
 id: string;
 status: "pending" | "running" | "complete" | "error";
 messages: BaseMessage[];
 result: string | undefined;
 toolCall: {
 id: string;
 name: string;
 args: {
 description: string;
 subagent_type: string;
 [key: string]: unknown;
 };
 };
 startedAt: number | undefined;
 completedAt: number | undefined;
}

```

| 
 | Property | Description
 | `id` | Unique identifier for this subagent instance
 | `status` | Lifecycle state: `pending` → `running` → `complete` or `error`
 | `messages` | The subagent’s own message stream, updated in real time
 | `result` | The final output text, available only when `status` is `"complete"`
 | `toolCall` | The tool call that spawned this subagent, including task metadata
 | `toolCall.args.description` | The task description the coordinator assigned to this subagent
 | `toolCall.args.subagent_type` | The type or name of the specialist (e.g., `"researcher"`, `"analyst"`)
 | `startedAt` | Timestamp when the subagent began executing
 | `completedAt` | Timestamp when the subagent finished

## [​](#linking-subagents-to-messages)Linking subagents to messages

The `getSubagentsByMessage` method returns the subagents spawned by a specific
AI message. This lets you render subagent cards directly beneath the
coordinator message that triggered them:

```
const turnSubagents = stream.getSubagentsByMessage(msg.id);

```

This returns an array of `SubagentStreamInterface` objects. If the message
didn’t spawn any subagents, it returns an empty array.

## [​](#building-the-subagentcard)Building the SubagentCard

Each subagent card shows the specialist’s name, task description, streaming
content or final result, and timing information:

```
import { AIMessage } from "@langchain/core/messages";

function SubagentCard({
 subagent,
}: {
 subagent: SubagentStreamInterface;
}) {
 const [expanded, setExpanded] = useState(true);

 const title =
 subagent.toolCall?.args?.subagent_type ?? `Agent ${subagent.id}`;
 const description = subagent.toolCall?.args?.description ?? "";

 const lastAIMessage = subagent.messages
 .filter(AIMessage.isInstance)
 .at(-1);

 const displayContent =
 subagent.status === "complete"
 ? subagent.result
 : typeof lastAIMessage?.content === "string"
 ? lastAIMessage.content
 : "";

 const elapsed = getElapsedTime(subagent.startedAt, subagent.completedAt);

 return (
 <div className="rounded-lg border bg-white shadow-sm">
 <button
 onClick={() => setExpanded(!expanded)}
 className="flex w-full items-center justify-between p-4"
 >
 <div className="flex items-center gap-3">
 <StatusIcon status={subagent.status} />
 <div>
 <h4 className="font-semibold capitalize">{title}</h4>
 <p className="text-xs text-gray-500">{description}</p>
 </div>
 </div>
 <div className="flex items-center gap-2">
 {elapsed && (
 <span className="text-xs text-gray-400">{elapsed}</span>
 )}
 <StatusBadge status={subagent.status} />
 </div>
 </button>

 {expanded && displayContent && (
 <div className="border-t px-4 py-3">
 <div className="prose prose-sm max-w-none line-clamp-6">
 {displayContent}
 {subagent.status === "running" && (
 <span className="inline-block h-4 w-1 animate-pulse bg-blue-500" />
 )}
 </div>
 </div>
 )}
 </div>
 );
}

function getElapsedTime(
 startedAt: number | undefined,
 completedAt: number | undefined
): string | null {
 if (!startedAt) return null;
 const end = completedAt ?? Date.now();
 const seconds = Math.round((end - startedAt) / 1000);
 if (seconds < 60) return `${seconds}s`;
 return `${Math.floor(seconds / 60)}m ${seconds % 60}s`;
}

```

## [​](#status-icons-and-badges)Status icons and badges

Consistent visual indicators help users parse subagent status at a glance:

```
function StatusIcon({ status }: { status: SubagentStreamInterface["status"] }) {
 switch (status) {
 case "pending":
 return <span className="text-gray-400">○</span>;
 case "running":
 return <span className="animate-spin text-blue-500">◉</span>;
 case "complete":
 return <span className="text-green-500">✓</span>;
 case "error":
 return <span className="text-red-500">✕</span>;
 }
}

function StatusBadge({ status }: { status: SubagentStreamInterface["status"] }) {
 const styles = {
 pending: "bg-gray-100 text-gray-600",
 running: "bg-blue-100 text-blue-700",
 complete: "bg-green-100 text-green-700",
 error: "bg-red-100 text-red-700",
 };

 return (
 <span className={`rounded-full px-2 py-0.5 text-xs font-medium ${styles[status]}`}>
 {status}
 </span>
 );
}

```

## [​](#progress-tracking)Progress tracking

Show a progress bar and counter so users know how many subagents have finished:

```
function SubagentProgress({
 subagents,
}: {
 subagents: SubagentStreamInterface[];
}) {
 const completed = subagents.filter((s) => s.status === "complete").length;
 const total = subagents.length;
 const percentage = total > 0 ? Math.round((completed / total) * 100) : 0;

 return (
 <div className="space-y-1">
 <div className="flex items-center justify-between text-xs text-gray-500">
 <span>Subagent progress</span>
 <span>
 {completed}/{total} complete
 </span>
 </div>
 <div className="h-2 overflow-hidden rounded-full bg-gray-200">
 <div
 className="h-full rounded-full bg-blue-500 transition-all duration-300"
 style={{ width: `${percentage}%` }}
 />
 </div>
 </div>
 );
}

```

## [​](#rendering-messages-with-subagent-cards)Rendering messages with subagent cards

The key layout pattern is to render each coordinator message, and if that message
spawned subagents, render their cards immediately below it:

```
function MessageWithSubagents({
 message,
 subagents,
}: {
 message: BaseMessage;
 subagents: SubagentStreamInterface[];
}) {
 if (message.type === "human") {
 return <HumanMessage content={message.content} />;
 }

 return (
 <div className="space-y-3">
 {message.content && (
 <div className="prose prose-sm max-w-none">
 {message.content}
 </div>
 )}

 {subagents.length > 0 && (
 <div className="ml-4 space-y-3 border-l-2 border-blue-200 pl-4">
 <SubagentProgress subagents={subagents} />
 {subagents.map((subagent) => (
 <SubagentCard key={subagent.id} subagent={subagent} />
 ))}
 </div>
 )}
 </div>
 );
}

```

## [​](#synthesis-indicator)Synthesis indicator

After all subagents complete, the coordinator takes time to synthesize their
results into a final response. Show a clear indicator during this phase:

```
function SynthesisIndicator({
 subagents,
 isLoading,
}: {
 subagents: SubagentStreamInterface[];
 isLoading: boolean;
}) {
 const allComplete =
 subagents.length > 0 &&
 subagents.every((s) => s.status === "complete" || s.status === "error");

 if (!allComplete || !isLoading) return null;

 return (
 <div className="flex items-center gap-2 rounded-lg bg-purple-50 px-4 py-2 text-sm text-purple-700">
 <span className="animate-spin">⟳</span>
 Synthesizing results from {subagents.length} subagent
 {subagents.length !== 1 ? "s" : ""}...
 </div>
 );
}

```

The synthesis phase can take several seconds for complex multi-expert
workflows. A clear “Synthesizing results…” indicator prevents users from
thinking the agent has stalled.

## [​](#debug-unfiltered-output)Debug unfiltered output

During development, you can temporarily set `filterSubagentMessages: false` to
see the raw, interleaved output from all subagents in the main message stream.
This is useful for verifying that subagent tokens are flowing correctly, but
should not be used in production UIs.

## [​](#use-cases)Use cases

Deep agent subagent cards are the right choice when your agent workflow
involves:

{' ' * (self.list_depth - 1)}- **Deep research** where a coordinator dispatches researchers to investigate
different facets of a question, then synthesizes their findings

{' ' * (self.list_depth - 1)}- **Multi-expert analysis** such as domain specialists (legal, financial, technical)
each contribute their perspective

{' ' * (self.list_depth - 1)}- **Complex task decomposition** where a planner breaks a large task into subtasks
and assigns each to a specialist worker

{' ' * (self.list_depth - 1)}- **Code review pipelines** where separate agents handle security review, style
checking, performance analysis, and documentation review

## [​](#accessing-the-full-subagents-map)Accessing the full subagents map

Beyond per-message lookup, you can access all subagents at once through
`stream.subagents`:

```
const allSubagents = [...stream.subagents.values()];
const running = allSubagents.filter((s) => s.status === "running");
const completed = allSubagents.filter((s) => s.status === "complete");
const errors = allSubagents.filter((s) => s.status === "error");

```

This is useful for building global progress indicators or dashboards that
summarize all subagent activity regardless of which coordinator message spawned
them.

## [​](#best-practices)Best practices

{' ' * (self.list_depth - 1)}- **Always set `filterSubagentMessages: true`**. Unfiltered streams produce an
unreadable interleaving of coordinator and subagent tokens.

{' ' * (self.list_depth - 1)}- **Show task descriptions**. The `toolCall.args.description` field tells users
exactly what each subagent was asked to do. Always display this prominently.

{' ' * (self.list_depth - 1)}- **Use collapsible cards**. In workflows with 5+ subagents, auto-collapse
completed cards so users can focus on active work.

{' ' * (self.list_depth - 1)}- **Display timing data**. Showing how long each subagent took helps users
understand performance characteristics and identify bottlenecks.

{' ' * (self.list_depth - 1)}- **Set an appropriate recursion limit**. Deep agent workflows with nested
subgraphs need higher limits than the default 25. Start with 100.

{' ' * (self.list_depth - 1)}- **Handle errors per subagent**. One subagent failing shouldn’t crash the
entire UI. Show the error in that subagent’s card while others continue
running.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/deepagents/frontend/subagent-streaming.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[OverviewPrevious](/oss/python/deepagents/frontend/overview)[Todo listNext](/oss/python/deepagents/frontend/todo-list)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
