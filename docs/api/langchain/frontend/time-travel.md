# Time travel

[Frontend](/oss/python/langchain/frontend/overview)[Patterns](/oss/python/langchain/frontend/markdown-messages)

# Time travel
Copy page

Inspect, navigate, and resume from any checkpoint in the conversation historyCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Every state change in a LangGraph agent creates a **checkpoint**, a complete
snapshot of the agent’s state at that moment. Time travel lets you inspect any
checkpoint, view the exact state the agent held, and **resume execution from
that point** to explore alternative paths. It’s a debugger, an undo button, and
an audit log all in one.

This feature requires the [LangGraph Agent Server](/langsmith/local-server). Run your agent locally with `langgraph dev` or [deploy it to LangSmith](/langsmith/deployment) to use this pattern.

## [​](#how-checkpoints-work)How checkpoints work

LangGraph persists agent state after every node execution. Each persisted state
is a `ThreadState` object that captures:

{' ' * (self.list_depth - 1)}- **checkpoint**: metadata identifying this specific snapshot (ID, timestamp)

{' ' * (self.list_depth - 1)}- **values**: the full agent state at this point (messages, custom keys)

{' ' * (self.list_depth - 1)}- **tasks**: the graph nodes that were scheduled to run next

{' ' * (self.list_depth - 1)}- **next**: the names of upcoming nodes in the execution plan

This creates a linear timeline of every decision the agent made, every tool it
called, and every response it produced. Your UI can render this timeline and let
users jump to any point.

## [​](#setting-up-usestream)Setting up useStream

Enable checkpoint history by passing `fetchStateHistory: true` to `useStream`.
This tells the hook to load the full checkpoint timeline for the current thread.
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

export function TimeTravelChat() {
 const stream = useStream<typeof myAgent>({
 apiUrl: AGENT_URL,
 assistantId: "time_travel",
 fetchStateHistory: true,
 });

 const history = stream.history ?? [];

 return (
 <div className="flex h-screen">
 <ChatPanel messages={stream.messages} />
 <TimelineSidebar
 history={history}
 onSelect={(cp) => stream.submit(null, { checkpoint: cp.checkpoint })}
 />
 </div>
 );
}

```

## [​](#the-threadstate-object)The ThreadState object

Each entry in the `history` array is a `ThreadState` representing one checkpoint
in the timeline:

```
interface ThreadState {
 checkpoint: {
 checkpoint_id: string;
 checkpoint_ns: string;
 };
 values: Record<string, unknown>;
 tasks: Array<{
 id: string;
 name: string;
 interrupts?: unknown[];
 }>;
 next: string[];
}

```

| 
 | Property | Description
 | `checkpoint` | Identifies this snapshot. Pass it to `submit` to resume from here
 | `values` | The complete agent state at this point, including `messages` and any custom state keys
 | `tasks` | The graph nodes that ran at this checkpoint, including their names and any interrupts
 | `next` | Names of nodes scheduled to execute after this checkpoint

## [​](#building-a-checkpoint-timeline)Building a checkpoint timeline

The timeline sidebar shows every checkpoint as a clickable entry. Each entry
displays the node that ran and how many messages existed at that point:

```
function TimelineSidebar({
 history,
 onSelect,
}: {
 history: ThreadState[];
 onSelect: (cp: ThreadState) => void;
}) {
 return (
 <aside className="w-80 overflow-y-auto border-l bg-gray-50 p-4">
 <h2 className="mb-4 text-sm font-semibold uppercase text-gray-500">
 Checkpoint Timeline
 </h2>
 <div className="space-y-2">
 {history.map((cp, i) => {
 const taskName = cp.tasks?.[0]?.name ?? "unknown";
 const msgCount = (cp.values?.messages as unknown[])?.length ?? 0;

 return (
 <button
 key={cp.checkpoint.checkpoint_id}
 onClick={() => onSelect(cp)}
 className="w-full rounded-lg border bg-white p-3 text-left
 hover:border-blue-400 hover:shadow-sm transition-all"
 >
 <div className="flex items-center justify-between">
 <span className="text-xs text-gray-400">#{i + 1}</span>
 <NodeBadge name={taskName} />
 </div>
 <p className="mt-1 text-sm font-medium">{taskName}</p>
 <p className="text-xs text-gray-500">
 {msgCount} message{msgCount !== 1 ? "s" : ""}
 </p>
 </button>
 );
 })}
 </div>
 </aside>
 );
}

```

## [​](#inspecting-checkpoint-state)Inspecting checkpoint state

Clicking a checkpoint should show the full state at that point. A JSON viewer
gives developers complete visibility into what the agent knew and decided:

```
function CheckpointInspector({ checkpoint }: { checkpoint: ThreadState }) {
 const [expanded, setExpanded] = useState(false);

 return (
 <div className="rounded-lg border bg-white p-4">
 <div className="flex items-center justify-between">
 <h3 className="font-semibold">
 Checkpoint {checkpoint.checkpoint.checkpoint_id.slice(0, 8)}...
 </h3>
 <button
 onClick={() => setExpanded(!expanded)}
 className="text-sm text-blue-600 hover:underline"
 >
 {expanded ? "Collapse" : "Expand"} state
 </button>
 </div>

 <div className="mt-2 space-y-1 text-sm">
 <p>
 <strong>Node:</strong>{" "}
 {checkpoint.tasks?.[0]?.name ?? "—"}
 </p>
 <p>
 <strong>Next:</strong>{" "}
 {checkpoint.next?.join(", ") || "—"}
 </p>
 <p>
 <strong>Messages:</strong>{" "}
 {(checkpoint.values?.messages as unknown[])?.length ?? 0}
 </p>
 </div>

 {expanded && (
 <div className="mt-3 max-h-96 overflow-auto rounded bg-gray-900 p-3">
 <pre className="text-xs text-gray-200">
 {JSON.stringify(checkpoint.values, null, 2)}
 </pre>
 </div>
 )}
 </div>
 );
}

```

For production UIs, consider using a proper JSON viewer component with
collapsible nodes instead of raw `JSON.stringify`. Libraries like
`react-json-view` or `react-json-tree` give users a much better exploration
experience.

## [​](#resuming-from-a-checkpoint)Resuming from a checkpoint

The core of time travel is the ability to **resume execution from any prior
checkpoint**. When a user selects a checkpoint, call `submit` with `null` input
and pass the checkpoint reference:

```
stream.submit(null, { checkpoint: selectedCheckpoint.checkpoint });

```

This tells LangGraph to:

{' ' * (self.list_depth - 1)}- Roll back to the selected checkpoint’s state

{' ' * (self.list_depth - 1)}- Re-execute the graph from that point forward

{' ' * (self.list_depth - 1)}- Stream the new results to the client

The existing messages after the selected checkpoint are replaced by the new
execution path. This effectively creates a **branch** in the conversation
timeline.
Resuming from a checkpoint does not delete the original timeline. The previous
checkpoints remain available in the history. This means users can always go back
and try a different path without losing any prior work.

## [​](#the-splitview-layout)The SplitView layout

Time travel works best with a split layout, with the main chat on the left and the
timeline on the right:

```
function TimeTravelLayout() {
 const stream = useStream<typeof myAgent>({
 apiUrl: AGENT_URL,
 assistantId: "time_travel",
 fetchStateHistory: true,
 });

 const [selectedCheckpoint, setSelectedCheckpoint] =
 useState<ThreadState | null>(null);

 const history = stream.history ?? [];

 return (
 <div className="flex h-screen">
 {/* Main chat area */}
 <main className="flex-1 overflow-y-auto p-6">
 <div className="mx-auto max-w-2xl space-y-4">
 {stream.messages.map((msg) => (
 <Message key={msg.id} message={msg} />
 ))}
 </div>
 <ChatInput
 onSubmit={(text) =>
 stream.submit({ messages: [{ type: "human", content: text }] })
 }
 isLoading={stream.isLoading}
 />
 </main>

 {/* Timeline sidebar */}
 <aside className="w-96 overflow-y-auto border-l bg-gray-50">
 <TimelineSidebar
 history={history}
 selected={selectedCheckpoint}
 onSelect={setSelectedCheckpoint}
 onResume={(cp) =>
 stream.submit(null, { checkpoint: cp.checkpoint })
 }
 />
 {selectedCheckpoint && (
 <CheckpointInspector checkpoint={selectedCheckpoint} />
 )}
 </aside>
 </div>
 );
}

```

## [​](#extracting-checkpoint-metadata)Extracting checkpoint metadata

Transform raw checkpoint data into display-friendly entries for your timeline:

```
function formatCheckpoints(history: ThreadState[]) {
 return history.map((cp, index) => ({
 index,
 id: cp.checkpoint?.checkpoint_id,
 taskName: cp.tasks?.[0]?.name ?? "unknown",
 messageCount: (cp.values?.messages as unknown[])?.length ?? 0,
 hasInterrupts: cp.tasks?.some((t) => t.interrupts?.length) ?? false,
 nextNodes: cp.next ?? [],
 }));
}

```

This makes it easy to render timeline entries with meaningful labels instead of
raw IDs.

## [​](#use-cases)Use cases

Time travel is invaluable across many scenarios:

{' ' * (self.list_depth - 1)}- **Debugging agent behavior**: step through the agent’s decisions to
understand why it chose a particular path

{' ' * (self.list_depth - 1)}- **Undoing actions**: if the agent took a wrong turn, resume from an earlier
checkpoint and try again

{' ' * (self.list_depth - 1)}- **Exploring alternatives**: fork from a mid-conversation checkpoint to see
how different inputs change the outcome

{' ' * (self.list_depth - 1)}- **Auditing**: review the complete history of an agent’s actions for
compliance, quality assurance, or post-incident analysis

{' ' * (self.list_depth - 1)}- **Teaching**: walk through an agent’s execution step by step to explain how
multi-step reasoning works

Time travel is especially powerful when combined with
[human-in-the-loop](/oss/python/langchain/frontend/human-in-the-loop) patterns. If a human reviewer
rejects an agent’s action at an interrupt, they can resume from the checkpoint
before the action was taken and provide corrective input.

## [​](#handling-interrupts-in-the-timeline)Handling interrupts in the timeline

Checkpoints that contain interrupts (human-in-the-loop pauses) deserve special
visual treatment. They represent moments where the agent stopped and waited for
human input:

```
function TimelineEntry({
 checkpoint,
 index,
}: {
 checkpoint: ThreadState;
 index: number;
}) {
 const hasInterrupt = checkpoint.tasks?.some(
 (t) => t.interrupts && t.interrupts.length > 0
 );

 return (
 <div
 className={`rounded-lg border p-3 ${
 hasInterrupt
 ? "border-amber-300 bg-amber-50"
 : "border-gray-200 bg-white"
 }`}
 >
 <div className="flex items-center gap-2">
 <span className="text-xs text-gray-400">#{index + 1}</span>
 {hasInterrupt && (
 <span className="rounded bg-amber-200 px-1.5 py-0.5 text-xs font-medium text-amber-800">
 Interrupt
 </span>
 )}
 </div>
 <p className="mt-1 text-sm font-medium">
 {checkpoint.tasks?.[0]?.name ?? "—"}
 </p>
 </div>
 );
}

```

## [​](#best-practices)Best practices

{' ' * (self.list_depth - 1)}- **Load history lazily**: for threads with hundreds of checkpoints, paginate
or load only the most recent N entries to keep the UI responsive.

{' ' * (self.list_depth - 1)}- **Show meaningful labels**: display node names and message counts instead of
raw checkpoint IDs. Users need context, not UUIDs.

{' ' * (self.list_depth - 1)}- **Confirm before resuming**: resuming from an old checkpoint replaces the
current execution path. Show a confirmation dialog so users don’t
accidentally lose the current conversation state.

{' ' * (self.list_depth - 1)}- **Highlight the current checkpoint**: make it visually obvious which
checkpoint corresponds to the current state of the conversation.

{' ' * (self.list_depth - 1)}- **Support keyboard navigation**: power users will want to step through
checkpoints with arrow keys. Add keyboard handlers to the timeline for a
smooth debugging experience.

{' ' * (self.list_depth - 1)}- **Diff state between checkpoints**: for advanced users, showing what changed
between two consecutive checkpoints can reveal exactly how the agent’s state
evolved at each step.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/frontend/time-travel.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Join & rejoin streamsPrevious](/oss/python/langchain/frontend/join-rejoin)[Generative UINext](/oss/python/langchain/frontend/generative-ui)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
