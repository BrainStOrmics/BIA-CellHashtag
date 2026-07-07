# Join &amp; rejoin streams

[Frontend](/oss/python/langchain/frontend/overview)[Patterns](/oss/python/langchain/frontend/markdown-messages)

# Join & rejoin streams
Copy page

Disconnect from and reconnect to running agent streamsCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Join and rejoin lets you disconnect from a running agent stream without stopping the agent, then reconnect to it later. The agent continues executing server-side while the client is away, and you pick up the stream exactly where you left off.

This feature requires the [LangGraph Agent Server](/langsmith/local-server). Run your agent locally with `langgraph dev` or [deploy it to LangSmith](/langsmith/deployment) to use this pattern.

## [​](#why-join-&-rejoin)Why join & rejoin?

Traditional streaming APIs tightly couple the client and server: if the client disconnects, the stream is lost. Join and rejoin breaks this coupling, enabling several important patterns:

{' ' * (self.list_depth - 1)}- **Network interruptions**: mobile users moving between cell towers or Wi-Fi networks can seamlessly resume

{' ' * (self.list_depth - 1)}- **Page navigation**: users navigating away from a chat page and returning later without losing progress

{' ' * (self.list_depth - 1)}- **Mobile backgrounding**: apps suspended by the OS can rejoin the stream when foregrounded

{' ' * (self.list_depth - 1)}- **Long-running tasks**: agents performing multi-minute operations (research, code generation, data analysis) where users don’t need to keep the page open

{' ' * (self.list_depth - 1)}- **Multi-device handoff**: start a conversation on your phone, rejoin on your desktop

## [​](#core-concepts)Core concepts

The join/rejoin pattern involves three key mechanisms:

| 
 | Method / Option | Purpose
 | `stream.stop()` | Disconnect the client from the stream without stopping the agent
 | `stream.joinStream(runId)` | Reconnect to an existing stream by its run ID
 | `onDisconnect: "continue"` | Submit option that tells the server to keep running after client disconnects
 | `streamResumable: true` | Submit option that enables the stream to be rejoined later
`stream.stop()` is fundamentally different from cancelling a run. Stopping only disconnects the **client**. The agent continues processing server-side. To actually cancel the agent’s execution, you would use interrupt or cancel mechanisms instead.

## [​](#setting-up-usestream)Setting up `useStream`

The key setup step is capturing the `run_id` from the `onCreated` callback so you can rejoin later.
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
import { useState } from "react";

function Chat() {
 const [savedRunId, setSavedRunId] = useState<string | null>(null);

 const stream = useStream<typeof myAgent>({
 apiUrl: "http://localhost:2024",
 assistantId: "join_rejoin",
 onCreated(run) {
 setSavedRunId(run.run_id);
 },
 });

 const isConnected = stream.isLoading;

 return (
 <div>
 <ConnectionStatus connected={isConnected} />
 <MessageList messages={stream.messages} />
 <ChatControls
 stream={stream}
 savedRunId={savedRunId}
 isConnected={isConnected}
 />
 </div>
 );
}

```

## [​](#submitting-with-resumable-options)Submitting with resumable options

When you submit a message, pass `onDisconnect: "continue"` and `streamResumable: true` to enable the join/rejoin flow:

```
stream.submit(
 { messages: [{ type: "human", content: text }] },
 {
 onDisconnect: "continue",
 streamResumable: true,
 }
);

```

| 
 | Option | Default | Description
 | `onDisconnect` | `"cancel"` | What happens when the client disconnects. `"continue"` keeps the agent running; `"cancel"` stops it.
 | `streamResumable` | `false` | When `true`, the server retains the stream state so a client can rejoin later.
Always use both options together. Setting `onDisconnect: "continue"` without `streamResumable: true` means the agent keeps running but you cannot rejoin the stream to see its output.

## [​](#disconnecting-from-a-stream)Disconnecting from a stream

Call `stream.stop()` to disconnect the client. The agent continues processing server-side.

```
stream.stop();

```

After calling `stop()`:

{' ' * (self.list_depth - 1)}- `stream.isLoading` becomes `false`

{' ' * (self.list_depth - 1)}- The message list retains all messages received up to the disconnect point

{' ' * (self.list_depth - 1)}- The agent continues running on the server

{' ' * (self.list_depth - 1)}- No new messages are received until you rejoin

## [​](#rejoining-a-stream)Rejoining a stream

Call `stream.joinStream(runId)` with the saved run ID to reconnect:

```
stream.joinStream(savedRunId);

```

After rejoining:

{' ' * (self.list_depth - 1)}- `stream.isLoading` becomes `true` again

{' ' * (self.list_depth - 1)}- Any messages generated while disconnected are delivered

{' ' * (self.list_depth - 1)}- New streaming messages resume in real-time

{' ' * (self.list_depth - 1)}- If the agent has already finished, you receive the final state immediately

## [​](#building-a-connection-status-indicator)Building a connection status indicator

A visual indicator helps users understand whether they are actively receiving updates from the agent.

```
function ConnectionStatus({ connected }: { connected: boolean }) {
 return (
 <div className="connection-status">
 <span
 className={`status-dot ${connected ? "connected" : "disconnected"}`}
 />
 <span className="status-text">
 {connected ? "Connected" : "Disconnected"}
 </span>
 </div>
 );
}

```

Style the indicator with a green/red dot:

```
.status-dot {
 width: 8px;
 height: 8px;
 border-radius: 50%;
 display: inline-block;
 margin-right: 6px;
}

.status-dot.connected {
 background-color: #22c55e;
 box-shadow: 0 0 4px #22c55e;
}

.status-dot.disconnected {
 background-color: #ef4444;
 box-shadow: 0 0 4px #ef4444;
}

```

## [​](#disconnect-and-rejoin-controls)Disconnect and rejoin controls

Provide explicit buttons for disconnecting and rejoining so users have full control:

```
function ChatControls({ stream, savedRunId, isConnected }) {
 const [input, setInput] = useState("");

 const handleSend = () => {
 if (!input.trim()) return;
 stream.submit(
 { messages: [{ type: "human", content: input.trim() }] },
 { onDisconnect: "continue", streamResumable: true }
 );
 setInput("");
 };

 return (
 <div className="controls">
 <div className="input-row">
 <input
 value={input}
 onChange={(e) => setInput(e.target.value)}
 placeholder="Type a message..."
 onKeyDown={(e) => e.key === "Enter" && handleSend()}
 />
 <button onClick={handleSend}>Send</button>
 </div>

 <div className="stream-controls">
 {isConnected ? (
 <button onClick={() => stream.stop()} className="disconnect-btn">
 Disconnect
 </button>
 ) : (
 savedRunId && (
 <button
 onClick={() => stream.joinStream(savedRunId)}
 className="rejoin-btn"
 >
 Rejoin stream
 </button>
 )
 )}
 </div>
 </div>
 );
}

```

## [​](#persisting-the-run-id)Persisting the run ID

For cross-session rejoin (e.g., the user closes the browser and returns later), persist the run ID to storage:

```
const stream = useStream<typeof myAgent>({
 apiUrl: "http://localhost:2024",
 assistantId: "join_rejoin",
 onCreated(run) {
 localStorage.setItem("activeRunId", run.run_id);
 },
});

// On page load, check for an active run
const existingRunId = localStorage.getItem("activeRunId");
if (existingRunId) {
 stream.joinStream(existingRunId);
}

```

Persisted run IDs should be cleaned up when a run completes. Listen for the stream to finish and remove the stored ID to avoid attempting to rejoin completed runs.

## [​](#error-handling)Error handling

Rejoining can fail if the run has expired, been deleted, or if the server has restarted. Handle these cases gracefully:

```
try {
 stream.joinStream(savedRunId);
} catch (error) {
 console.error("Failed to rejoin stream:", error);
 // Clear stale run ID and inform the user
 setSavedRunId(null);
 localStorage.removeItem("activeRunId");
}

```

## [​](#complete-example)Complete example

```
function JoinRejoinChat() {
 const [savedRunId, setSavedRunId] = useState<string | null>(null);
 const [input, setInput] = useState("");

 const stream = useStream<typeof myAgent>({
 apiUrl: "http://localhost:2024",
 assistantId: "join_rejoin",
 onCreated(run) {
 setSavedRunId(run.run_id);
 },
 });

 const isConnected = stream.isLoading;

 const handleSend = () => {
 if (!input.trim()) return;
 stream.submit(
 { messages: [{ type: "human", content: input.trim() }] },
 { onDisconnect: "continue", streamResumable: true }
 );
 setInput("");
 };

 return (
 <div className="chat-container">
 <header>
 <h2>Join & Rejoin Demo</h2>
 <ConnectionStatus connected={isConnected} />
 </header>

 <div className="messages">
 {stream.messages.map((msg, i) => (
 <MessageBubble key={i} message={msg} />
 ))}
 </div>

 <div className="controls">
 <form onSubmit={(e) => { e.preventDefault(); handleSend(); }}>
 <input
 value={input}
 onChange={(e) => setInput(e.target.value)}
 placeholder="Type a message..."
 />
 <button type="submit">Send</button>
 </form>

 <div className="stream-actions">
 {isConnected ? (
 <button onClick={() => stream.stop()}>
 Disconnect
 </button>
 ) : (
 savedRunId && (
 <button onClick={() => stream.joinStream(savedRunId)}>
 Rejoin stream
 </button>
 )
 )}
 </div>
 </div>
 </div>
 );
}

```

## [​](#best-practices)Best practices

{' ' * (self.list_depth - 1)}- **Always save the run ID**: without it, rejoining is impossible. Use both component state and persistent storage for resilience.

{' ' * (self.list_depth - 1)}- **Show clear connection state**: users should always know whether they are receiving live updates or viewing a snapshot.

{' ' * (self.list_depth - 1)}- **Auto-rejoin on visibility change**: use the Page Visibility API to automatically rejoin when the user returns to the tab.

{' ' * (self.list_depth - 1)}- **Set reasonable timeouts**: if a rejoin attempt takes too long, fall back to fetching the thread history instead.

{' ' * (self.list_depth - 1)}- **Clean up completed runs**: remove persisted run IDs when the agent finishes to avoid stale rejoin attempts.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/frontend/join-rejoin.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Message queuesPrevious](/oss/python/langchain/frontend/message-queues)[Time travelNext](/oss/python/langchain/frontend/time-travel)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
